#include "reconstruction.h"

#include <algorithm>
#include <cmath>

#include "config.h"
#include "kernels.h"
#include "math_utils.h"
#include "mfm.h"

namespace Reconstruction {

namespace {  // Anonymous namespace for private helper functions

inline double compute_gradient_alpha(double d_phi_max, double d_phi_min,
                                     double d_mid_max, double d_mid_min,
                                     double beta) {
    double alpha = 1.0;
    if (d_mid_max > 0.0) {
        if (d_phi_max <= 0.0)
            alpha = 0.0;
        else
            alpha = std::min(alpha, beta * d_phi_max / d_mid_max);
    }
    if (d_mid_min < 0.0) {
        if (d_phi_min >= 0.0)
            alpha = 0.0;
        else
            alpha = std::min(alpha, beta * d_phi_min / d_mid_min);
    }
    return alpha;
}

inline int math_sign(double x) { return (x > 0.0) ? 1 : ((x < 0.0) ? -1 : 0); }

double apply_pairwise_limiter(double phi_L_center, double phi_R_center,
                              double phi_mid_0, double phi_bar,
                              const Config& config) {
    double d_phi = std::abs(phi_L_center - phi_R_center);
    if (d_phi < 1e-14) return phi_L_center;

    double phi_min = std::min(phi_L_center, phi_R_center);
    double phi_max = std::max(phi_L_center, phi_R_center);

    double delta_1 = config.pairwise_limiter_psi1 * d_phi;
    double delta_2 = config.pairwise_limiter_psi2 * d_phi;

    double phi_minus =
        (math_sign(phi_min - delta_1) == math_sign(phi_min) || phi_min == 0.0)
            ? phi_min - delta_1
            : phi_min / (1.0 + delta_1 / std::abs(phi_min));
    double phi_plus =
        (math_sign(phi_max + delta_1) == math_sign(phi_max) || phi_max == 0.0)
            ? phi_max + delta_1
            : phi_max / (1.0 + delta_1 / std::abs(phi_max));

    if (phi_L_center < phi_R_center) {
        return std::max(phi_minus, std::min(phi_bar + delta_2, phi_mid_0));
    } else {
        return std::min(phi_plus, std::max(phi_bar - delta_2, phi_mid_0));
    }
}

inline double compute_velocity_alpha(double vi, double vj, double dv_raw,
                                     double fraction, const Config& config) {
    if (std::abs(dv_raw) < 1e-14) return 1.0;

    double psi_1 = config.pairwise_limiter_psi1;
    double psi_2 = config.pairwise_limiter_psi2;

    double d_v = std::abs(vi - vj);
    double v_min = std::min(vi, vj) - psi_1 * d_v;
    double v_max = std::max(vi, vj) + psi_1 * d_v;
    double v_bar = vi + fraction * (vj - vi);
    double v_L_raw = vi + dv_raw;
    double v_lim;

    if (vi < vj) {
        v_lim = std::max(v_min, std::min(v_bar + psi_2 * d_v, v_L_raw));
    } else {
        v_lim = std::min(v_max, std::max(v_bar - psi_2 * d_v, v_L_raw));
    }
    return std::max(0.0, std::min(1.0, (v_lim - vi) / dv_raw));
}

Eigen::Vector3d compute_face_area_vector(
    const ParticleState& p_i, bool ill_conditioned_i, const ParticleState& p_j,
    bool ill_conditioned_j, const Eigen::Matrix3d& B_matrix_i,
    const Eigen::Matrix3d& B_matrix_j, double r, const Eigen::Vector3d& dx_vec,
    bool& used_sph_fallback) {
    used_sph_fallback = false;

    double V_i = p_i.mass / p_i.rho;
    double V_j = p_j.mass / p_j.rho;

    double W_i, W_j;
    Kernels::cubic_spline_value(r, p_i.h, W_i);
    Kernels::cubic_spline_value(r, p_j.h, W_j);

#ifdef USE_GEOMETRIC_MATRIX
    Eigen::Vector3d Area_vec = (V_j * W_i * (B_matrix_i * dx_vec)) +
                               (V_i * W_j * (B_matrix_j * dx_vec));
#else
    Eigen::Vector3d Area_vec = (V_i * V_j * W_i * (B_matrix_i * dx_vec)) +
                               (V_i * V_j * W_j * (B_matrix_j * dx_vec));
#endif

    double facenormal_dot_dp = Area_vec.dot(dx_vec);

    if (facenormal_dot_dp < 0.0 || ill_conditioned_i || ill_conditioned_j) {
        used_sph_fallback = true;
        double dW_dr_i, dW_dr_j, dummy_phi;
        Kernels::adaptive_gravity_terms(r, p_i.h, dummy_phi, dW_dr_i);
        Kernels::adaptive_gravity_terms(r, p_j.h, dummy_phi, dW_dr_j);
        double face_area_mag = -(V_i * V_i * dW_dr_i + V_j * V_j * dW_dr_j) / r;
        return face_area_mag * dx_vec;
    }

    return Area_vec;
}

}  // end anonymous namespace

ParticleGradients compute_single_particle_gradients(
    size_t i, const int* neighbor_indices, size_t num_neighbors,
    const FluidStateArrays& data, const Eigen::Matrix3d& B,
    bool ill_conditioned, const Config& config) {
    const Vec3& domain_size = config.domain_size;

    ParticleGradients out;
    out.grad_rho = Eigen::Vector3d::Zero();
    out.grad_p = Eigen::Vector3d::Zero();
    out.grad_vx = Eigen::Vector3d::Zero();
    out.grad_vy = Eigen::Vector3d::Zero();
    out.grad_vz = Eigen::Vector3d::Zero();
    out.ill_conditioned = ill_conditioned;

    if (!ill_conditioned) {
        Eigen::Vector3d sum_rho = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_p = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vx = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vy = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vz = Eigen::Vector3d::Zero();

        double d_rho_max = 0.0, d_rho_min = 0.0;
        double d_p_max = 0.0, d_p_min = 0.0;
        double d_v_max = 0.0;

        for (size_t n = 0; n < num_neighbors; ++n) {
            int j = neighbor_indices[n];

            double dx = periodic_displacement(data.pos_x[j] - data.pos_x[i],
                                              domain_size.x);
            double dy = periodic_displacement(data.pos_y[j] - data.pos_y[i],
                                              domain_size.y);
            double dz = periodic_displacement(data.pos_z[j] - data.pos_z[i],
                                              domain_size.z);

            double r2 = dx * dx + dy * dy + dz * dz;
            double r = std::sqrt(r2);

            double dv_x = data.vel_x[j] - data.vel_x[i];
            double dv_y = data.vel_y[j] - data.vel_y[i];
            double dv_z = data.vel_z[j] - data.vel_z[i];
            double dv_mag = std::sqrt(dv_x * dv_x + dv_y * dv_y + dv_z * dv_z);

            // Track min/max bounds across the entire interaction list
            d_rho_max = std::max(d_rho_max, data.rho[j] - data.rho[i]);
            d_rho_min = std::min(d_rho_min, data.rho[j] - data.rho[i]);
            d_p_max = std::max(d_p_max, data.pressure[j] - data.pressure[i]);
            d_p_min = std::min(d_p_min, data.pressure[j] - data.pressure[i]);
            d_v_max = std::max(d_v_max, dv_mag);

            double W;
            Kernels::cubic_spline_value(r, data.h[i], W);
            Eigen::Vector3d dx_vec(dx, dy, dz);
            double V_j = data.mass[j] / data.rho[j];
#ifdef USE_GEOMETRIC_MATRIX
            double weight = W;
#else
            double weight = W * V_j;
#endif

            sum_rho += (data.rho[j] - data.rho[i]) * dx_vec * weight;
            sum_p += (data.pressure[j] - data.pressure[i]) * dx_vec * weight;
            sum_vx += dv_x * dx_vec * weight;
            sum_vy += dv_y * dx_vec * weight;
            sum_vz += dv_z * dx_vec * weight;
        }

        out.grad_rho = B * sum_rho;
        out.grad_p = B * sum_p;
        out.grad_vx = B * sum_vx;
        out.grad_vy = B * sum_vy;
        out.grad_vz = B * sum_vz;

        double phi_mid_max_rho = 0.0, phi_mid_min_rho = 0.0;
        double phi_mid_max_p = 0.0, phi_mid_min_p = 0.0;
        double phi_mid_max_v = 0.0;

        for (size_t n = 0; n < num_neighbors; ++n) {
            int j = neighbor_indices[n];

            double dx = periodic_displacement(data.pos_x[j] - data.pos_x[i],
                                              domain_size.x);
            double dy = periodic_displacement(data.pos_y[j] - data.pos_y[i],
                                              domain_size.y);
            double dz = periodic_displacement(data.pos_z[j] - data.pos_z[i],
                                              domain_size.z);

#ifdef USE_MIDPOINT_QUADRATURE
            double fraction_i = 0.5;
#else
            double fraction_i = data.h[i] / (data.h[i] + data.h[j]);
#endif
            Eigen::Vector3d dx_face_i =
                fraction_i * Eigen::Vector3d(dx, dy, dz);

            double d_mid_rho = out.grad_rho.dot(dx_face_i);
            phi_mid_max_rho = std::max(phi_mid_max_rho, d_mid_rho);
            phi_mid_min_rho = std::min(phi_mid_min_rho, d_mid_rho);

            double d_mid_p = out.grad_p.dot(dx_face_i);
            phi_mid_max_p = std::max(phi_mid_max_p, d_mid_p);
            phi_mid_min_p = std::min(phi_mid_min_p, d_mid_p);

            // Compute the extrapolated velocity difference vector at the face
            Eigen::Vector3d d_mid_v(out.grad_vx.dot(dx_face_i),
                                    out.grad_vy.dot(dx_face_i),
                                    out.grad_vz.dot(dx_face_i));

            phi_mid_max_v = std::max(phi_mid_max_v, d_mid_v.norm());
        }

        double beta = config.gradient_limiter_beta_max;
        if (data.cond_num[i] > 0.0) {
            beta = std::max(config.gradient_limiter_beta_min,
                            config.gradient_limiter_beta_max *
                                std::min(1.0, N_cond_crit / data.cond_num[i]));
        }

        out.grad_rho *= compute_gradient_alpha(
            d_rho_max, d_rho_min, phi_mid_max_rho, phi_mid_min_rho, beta);
        out.grad_p *= compute_gradient_alpha(d_p_max, d_p_min, phi_mid_max_p,
                                             phi_mid_min_p, beta);

        double alpha_v =
            compute_gradient_alpha(d_v_max, 0.f, phi_mid_max_v, 0.f, beta);

        out.grad_vx *= alpha_v;
        out.grad_vy *= alpha_v;
        out.grad_vz *= alpha_v;

    } else {
        out.grad_rho.setZero();
        out.grad_p.setZero();
        out.grad_vx.setZero();
        out.grad_vy.setZero();
        out.grad_vz.setZero();

        /*for (size_t n = 0; n < num_neighbors; ++n) {
            int j = neighbor_indices[n];

            double dx = periodic_displacement(data.pos_x[j] - data.pos_x[i],
                                              domain_size.x);
            double dy = periodic_displacement(data.pos_y[j] - data.pos_y[i],
                                              domain_size.y);
            double dz = periodic_displacement(data.pos_z[j] - data.pos_z[i],
                                              domain_size.z);
            double r2 = dx * dx + dy * dy + dz * dz;
            double r = std::sqrt(r2);
            double dphi_dr_dummy, dW_dr;
            Kernels::adaptive_gravity_terms(r, data.h[i], dphi_dr_dummy, dW_dr);

            double V_j = data.mass[j] / data.rho[j];
            Eigen::Vector3d dx_vec(dx, dy, dz);
            Eigen::Vector3d grad_W_i = -dW_dr * dx_vec / r;

            out.grad_rho += V_j * (data.rho[j] - data.rho[i]) * grad_W_i;
            out.grad_p +=
                V_j * (data.pressure[j] - data.pressure[i]) * grad_W_i;
            out.grad_vx += V_j * (data.vel_x[j] - data.vel_x[i]) * grad_W_i;
            out.grad_vy += V_j * (data.vel_y[j] - data.vel_y[i]) * grad_W_i;
            out.grad_vz += V_j * (data.vel_z[j] - data.vel_z[i]) * grad_W_i;
        }*/
    }

    return out;
}

ReconstructedFace compute_face_reconstruction(
    const ParticleState& p_i, const ParticleGradients& grad_i,
    const ParticleState& p_j, const ParticleGradients& grad_j,
    const Eigen::Matrix3d& B_matrix_i, const Eigen::Matrix3d& B_matrix_j,
    const Config& config, double density_floor, double pressure_floor) {
    const Vec3& domain_size = config.domain_size;
    ReconstructedFace face;
    face.is_valid = false;
    face.used_sph_fallback = false;

    double dx = periodic_displacement(p_j.pos.x() - p_i.pos.x(), domain_size.x);
    double dy = periodic_displacement(p_j.pos.y() - p_i.pos.y(), domain_size.y);
    double dz = periodic_displacement(p_j.pos.z() - p_i.pos.z(), domain_size.z);

    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < 1e-24) return face;

    face.r = std::sqrt(r2);
    Eigen::Vector3d dx_vec(dx, dy, dz);
    face.area_vec = compute_face_area_vector(
        p_i, grad_i.ill_conditioned, p_j, grad_j.ill_conditioned, B_matrix_i,
        B_matrix_j, face.r, dx_vec, face.used_sph_fallback);

    double A_mag = face.area_vec.norm();
    if (A_mag < 1e-20) return face;

    face.n = face.area_vec / A_mag;

    // Zeroth order
    if (face.used_sph_fallback) {
        face.rho_L = std::max(p_i.rho, density_floor);
        face.rho_R = std::max(p_j.rho, density_floor);

        face.p_L = std::max(p_i.pressure, pressure_floor);
        face.p_R = std::max(p_j.pressure, pressure_floor);

        face.v_L = p_i.vel;
        face.v_R = p_j.vel;

        face.is_valid = true;
        return face;
    }

#ifdef USE_MIDPOINT_QUADRATURE
    double fraction_i = 0.5;
#else
    double fraction_i = p_i.h / (p_i.h + p_j.h);
#endif
    double fraction_j = 1.0 - fraction_i;

    Eigen::Vector3d dx_face_i = fraction_i * dx_vec;
    Eigen::Vector3d dx_face_j = -fraction_j * dx_vec;

    double rho_bar = p_i.rho + fraction_i * (p_j.rho - p_i.rho);
    double p_bar = p_i.pressure + fraction_i * (p_j.pressure - p_i.pressure);

    face.rho_L = apply_pairwise_limiter(
        p_i.rho, p_j.rho, p_i.rho + grad_i.grad_rho.dot(dx_face_i), rho_bar,
        config);
    face.rho_R = apply_pairwise_limiter(
        p_j.rho, p_i.rho, p_j.rho + grad_j.grad_rho.dot(dx_face_j), rho_bar,
        config);

    face.p_L = apply_pairwise_limiter(
        p_i.pressure, p_j.pressure, p_i.pressure + grad_i.grad_p.dot(dx_face_i),
        p_bar, config);
    face.p_R = apply_pairwise_limiter(
        p_j.pressure, p_i.pressure, p_j.pressure + grad_j.grad_p.dot(dx_face_j),
        p_bar, config);

    // Galilean & Rotationally Invariant Velocity Limiting
    // Project center velocities onto the face normal
    double vn_i = p_i.vel.dot(face.n);
    double vn_j = p_j.vel.dot(face.n);

    // Left face velocity reconstruction
    Eigen::Vector3d dv_L_raw(grad_i.grad_vx.dot(dx_face_i),
                             grad_i.grad_vy.dot(dx_face_i),
                             grad_i.grad_vz.dot(dx_face_i));

    // Project the extrapolated velocity difference onto the normal
    double dvn_L_raw = dv_L_raw.dot(face.n);

    // Compute the rotationally invariant alpha
    double alpha_v_L =
        compute_velocity_alpha(vn_i, vn_j, dvn_L_raw, fraction_i, config);

    // Scale the full 3D vector gradient symmetrically
    face.v_L = p_i.vel + alpha_v_L * dv_L_raw;

    // Right face velocity reconstruction
    Eigen::Vector3d dv_R_raw(grad_j.grad_vx.dot(dx_face_j),
                             grad_j.grad_vy.dot(dx_face_j),
                             grad_j.grad_vz.dot(dx_face_j));

    double dvn_R_raw = dv_R_raw.dot(face.n);
    double alpha_v_R =
        compute_velocity_alpha(vn_j, vn_i, dvn_R_raw, fraction_j, config);

    face.v_R = p_j.vel + alpha_v_R * dv_R_raw;

    face.rho_L = std::max(face.rho_L, density_floor);
    face.rho_R = std::max(face.rho_R, density_floor);
    face.p_L = std::max(face.p_L, pressure_floor);
    face.p_R = std::max(face.p_R, pressure_floor);

    face.is_valid = true;
    return face;
}

}  // namespace Reconstruction
