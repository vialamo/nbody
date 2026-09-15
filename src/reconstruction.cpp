#include "reconstruction.h"

#include <algorithm>
#include <cmath>

#include "kernels.h"

// Enable to use the limiter explained in the Gizmo paper.
// Disable to use the limiter used in the Gizmo code.
// #define THEORETICAL_LIMITER

// Disable limiter completely
// #define DISABLE_LIMITER

#ifdef DISABLE_LIMITER
#define ZEROTH_ORDER_RECONSTRUCTION
#endif

namespace Reconstruction {

namespace {  // Anonymous namespace for private helper functions

inline double periodic_displacement(double dx, double domain_size) {
    if (dx > 0.5 * domain_size) return dx - domain_size;
    if (dx < -0.5 * domain_size) return dx + domain_size;
    return dx;
}

#ifdef THEORETICAL_LIMITER
inline int math_sign(double x) { return (x > 0.0) ? 1 : ((x < 0.0) ? -1 : 0); }

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

double apply_pairwise_limiter(double phi_L_center, double phi_R_center,
                              double phi_mid_0, double phi_bar) {
    double psi_1 = 0.5;
    double psi_2 = 0.25;

    double d_phi = std::abs(phi_L_center - phi_R_center);
    if (d_phi < 1e-14) return phi_L_center;

    double phi_min = std::min(phi_L_center, phi_R_center);
    double phi_max = std::max(phi_L_center, phi_R_center);

    double delta_1 = psi_1 * d_phi;
    double delta_2 = psi_2 * d_phi;

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
#else
inline void scalar_limiter(Eigen::Vector3d& grad, double valmax, double valmin,
                           double alim, double h, double shoot_tol,
                           bool pos_preserve, double d_max, double val_cen) {
    double d_abs = grad.norm();
    if (d_abs > 0.0) {
        double cfac = 1.0 / (alim * h * d_abs);
        double abs_max = std::abs(valmax);
        double abs_min = std::abs(valmin);
        if (abs_max < abs_min) std::swap(abs_max, abs_min);

        double f_corr_overshoot =
            std::min(abs_min + shoot_tol * abs_max, abs_max);
        cfac *= f_corr_overshoot;

        if (pos_preserve) {
            constexpr double MIN_REAL_NUMBER = 1e-30;
            double fmin = std::min(
                val_cen,
                std::max(0.0, std::max(MIN_REAL_NUMBER * val_cen,
                                       std::min(0.5 * (val_cen + valmin),
                                                val_cen - f_corr_overshoot))));
            cfac = std::min((((val_cen - fmin) / d_max) / d_abs), cfac);
        }
        if (cfac < 1.0) grad *= cfac;
    }
}
#endif

inline double compute_condition_number(const Eigen::Matrix3d& E,
                                       const Eigen::Matrix3d& B) {
    return (1.0 / 3.0) * std::sqrt(E.squaredNorm() * B.squaredNorm());
}

}  // end anonymous namespace

ParticleGradients compute_single_particle_gradients(
    const ParticleState& p_i, const std::vector<ParticleState>& neighbors,
    double domain_size) {
    ParticleGradients out;
    out.B_matrix = Eigen::Matrix3d::Zero();
    out.grad_rho = Eigen::Vector3d::Zero();
    out.grad_p = Eigen::Vector3d::Zero();
    out.grad_vx = Eigen::Vector3d::Zero();
    out.grad_vy = Eigen::Vector3d::Zero();
    out.grad_vz = Eigen::Vector3d::Zero();
    out.ill_conditioned = false;

    Eigen::Matrix3d E = Eigen::Matrix3d::Zero();

    // Build the E Matrix
    for (const auto& nj : neighbors) {
        double dx =
            periodic_displacement(nj.pos.x() - p_i.pos.x(), domain_size);
        double dy =
            periodic_displacement(nj.pos.y() - p_i.pos.y(), domain_size);
        double dz =
            periodic_displacement(nj.pos.z() - p_i.pos.z(), domain_size);

        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < p_i.h * p_i.h && r2 > 1e-24) {
            double r = std::sqrt(r2);
            double W;
            Kernels::cubic_spline_value(r, p_i.h, W);

            double V_j = nj.mass / nj.rho;
            double weight = W * V_j;

            E(0, 0) += dx * dx * weight;
            E(0, 1) += dx * dy * weight;
            E(0, 2) += dx * dz * weight;
            E(1, 1) += dy * dy * weight;
            E(1, 2) += dy * dz * weight;
            E(2, 2) += dz * dz * weight;
        }
    }

    E(1, 0) = E(0, 1);
    E(2, 0) = E(0, 2);
    E(2, 1) = E(1, 2);

    double det = E.determinant();
    out.ill_conditioned = true;
    out.condition_number = -1.0;

    if (std::abs(det) > 1e-30) {
        Eigen::Matrix3d temp_B = E.inverse();
        double N_cond = compute_condition_number(E, temp_B);
        out.condition_number = N_cond;

        constexpr double N_cond_crit = 1000.0;
        if (N_cond <= N_cond_crit) {
            out.B_matrix = temp_B;
            out.ill_conditioned = false;
        }
    }

    if (!out.ill_conditioned) {
        out.B_matrix = E.inverse();
        Eigen::Vector3d sum_rho = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_p = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vx = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vy = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vz = Eigen::Vector3d::Zero();

        double d_rho_max = 0.0, d_rho_min = 0.0;
        double d_p_max = 0.0, d_p_min = 0.0;
        double d_vx_max = 0.0, d_vx_min = 0.0;
        double d_vy_max = 0.0, d_vy_min = 0.0;
        double d_vz_max = 0.0, d_vz_min = 0.0;
        double r_max = 1e-12;

        for (const auto& nj : neighbors) {
            double dx =
                periodic_displacement(nj.pos.x() - p_i.pos.x(), domain_size);
            double dy =
                periodic_displacement(nj.pos.y() - p_i.pos.y(), domain_size);
            double dz =
                periodic_displacement(nj.pos.z() - p_i.pos.z(), domain_size);

            double r2 = dx * dx + dy * dy + dz * dz;
            if ((r2 < p_i.h * p_i.h || r2 < nj.h * nj.h) && r2 > 1e-24) {
                double r = std::sqrt(r2);
                r_max = std::max(r_max, r);

                d_rho_max = std::max(d_rho_max, nj.rho - p_i.rho);
                d_rho_min = std::min(d_rho_min, nj.rho - p_i.rho);
                d_p_max = std::max(d_p_max, nj.pressure - p_i.pressure);
                d_p_min = std::min(d_p_min, nj.pressure - p_i.pressure);
                d_vx_max = std::max(d_vx_max, nj.vel.x() - p_i.vel.x());
                d_vx_min = std::min(d_vx_min, nj.vel.x() - p_i.vel.x());
                d_vy_max = std::max(d_vy_max, nj.vel.y() - p_i.vel.y());
                d_vy_min = std::min(d_vy_min, nj.vel.y() - p_i.vel.y());
                d_vz_max = std::max(d_vz_max, nj.vel.z() - p_i.vel.z());
                d_vz_min = std::min(d_vz_min, nj.vel.z() - p_i.vel.z());
            }

            if (r2 < p_i.h * p_i.h && r2 > 1e-24) {
                double r = std::sqrt(r2);
                double W;
                Kernels::cubic_spline_value(r, p_i.h, W);

                double V_j = nj.mass / nj.rho;
                double weight = W * V_j;
                Eigen::Vector3d dx_vec(dx, dy, dz);

                sum_rho += (nj.rho - p_i.rho) * dx_vec * weight;
                sum_p += (nj.pressure - p_i.pressure) * dx_vec * weight;
                sum_vx += (nj.vel.x() - p_i.vel.x()) * dx_vec * weight;
                sum_vy += (nj.vel.y() - p_i.vel.y()) * dx_vec * weight;
                sum_vz += (nj.vel.z() - p_i.vel.z()) * dx_vec * weight;
            }
        }

        out.grad_rho = out.B_matrix * sum_rho;
        out.grad_p = out.B_matrix * sum_p;
        out.grad_vx = out.B_matrix * sum_vx;
        out.grad_vy = out.B_matrix * sum_vy;
        out.grad_vz = out.B_matrix * sum_vz;
        out.raw_sum_p = sum_p;

#ifndef DISABLE_LIMITER
#ifdef THEORETICAL_LIMITER
        double phi_mid_max_rho = 0.0, phi_mid_min_rho = 0.0;
        double phi_mid_max_p = 0.0, phi_mid_min_p = 0.0;
        double phi_mid_max_vx = 0.0, phi_mid_min_vx = 0.0;
        double phi_mid_max_vy = 0.0, phi_mid_min_vy = 0.0;
        double phi_mid_max_vz = 0.0, phi_mid_min_vz = 0.0;

        for (const auto& nj : neighbors) {
            double dx =
                periodic_displacement(nj.pos.x() - p_i.pos.x(), domain_size);
            double dy =
                periodic_displacement(nj.pos.y() - p_i.pos.y(), domain_size);
            double dz =
                periodic_displacement(nj.pos.z() - p_i.pos.z(), domain_size);

            double r2 = dx * dx + dy * dy + dz * dz;
            if ((r2 < p_i.h * p_i.h || r2 < nj.h * nj.h) && r2 > 1e-24) {
                double fraction_i = p_i.h / (p_i.h + nj.h);
                Eigen::Vector3d dx_face_i =
                    fraction_i * Eigen::Vector3d(dx, dy, dz);

                double d_mid_rho = out.grad_rho.dot(dx_face_i);
                phi_mid_max_rho = std::max(phi_mid_max_rho, d_mid_rho);
                phi_mid_min_rho = std::min(phi_mid_min_rho, d_mid_rho);

                double d_mid_p = out.grad_p.dot(dx_face_i);
                phi_mid_max_p = std::max(phi_mid_max_p, d_mid_p);
                phi_mid_min_p = std::min(phi_mid_min_p, d_mid_p);

                double d_mid_vx = out.grad_vx.dot(dx_face_i);
                phi_mid_max_vx = std::max(phi_mid_max_vx, d_mid_vx);
                phi_mid_min_vx = std::min(phi_mid_min_vx, d_mid_vx);

                double d_mid_vy = out.grad_vy.dot(dx_face_i);
                phi_mid_max_vy = std::max(phi_mid_max_vy, d_mid_vy);
                phi_mid_min_vy = std::min(phi_mid_min_vy, d_mid_vy);

                double d_mid_vz = out.grad_vz.dot(dx_face_i);
                phi_mid_max_vz = std::max(phi_mid_max_vz, d_mid_vz);
                phi_mid_min_vz = std::min(phi_mid_min_vz, d_mid_vz);
            }
        }

        double beta = 1.5;
        out.grad_rho *= compute_gradient_alpha(
            d_rho_max, d_rho_min, phi_mid_max_rho, phi_mid_min_rho, beta);
        out.grad_p *= compute_gradient_alpha(d_p_max, d_p_min, phi_mid_max_p,
                                             phi_mid_min_p, beta);
        out.grad_vx *= compute_gradient_alpha(
            d_vx_max, d_vx_min, phi_mid_max_vx, phi_mid_min_vx, beta);
        out.grad_vy *= compute_gradient_alpha(
            d_vy_max, d_vy_min, phi_mid_max_vy, phi_mid_min_vy, beta);
        out.grad_vz *= compute_gradient_alpha(
            d_vz_max, d_vz_min, phi_mid_max_vz, phi_mid_min_vz, beta);
#else
        double alim = 0.5;
        double h_lim = std::max(p_i.h, r_max);
        double stol = 0.1;

        scalar_limiter(out.grad_rho, d_rho_max, d_rho_min, alim, h_lim, 0.0,
                       true, h_lim, p_i.rho);
        scalar_limiter(out.grad_p, d_p_max, d_p_min, alim, h_lim, stol, true,
                       h_lim, p_i.pressure);
        scalar_limiter(out.grad_vx, d_vx_max, d_vx_min, alim, h_lim, stol,
                       false, h_lim, p_i.vel.x());
        scalar_limiter(out.grad_vy, d_vy_max, d_vy_min, alim, h_lim, stol,
                       false, h_lim, p_i.vel.y());
        scalar_limiter(out.grad_vz, d_vz_max, d_vz_min, alim, h_lim, stol,
                       false, h_lim, p_i.vel.z());
#endif
#endif
    } else {
        out.B_matrix = Eigen::Matrix3d::Identity();
        out.ill_conditioned = true;

        double trace_E = E(0, 0) + E(1, 1) + E(2, 2);
        if (trace_E > 1e-24) {
            out.B_matrix = (3.0 / trace_E) * Eigen::Matrix3d::Identity();
        } else {
            out.B_matrix =
                (1.0 / (p_i.h * p_i.h)) * Eigen::Matrix3d::Identity();
        }

        for (const auto& nj : neighbors) {
            double dx =
                periodic_displacement(nj.pos.x() - p_i.pos.x(), domain_size);
            double dy =
                periodic_displacement(nj.pos.y() - p_i.pos.y(), domain_size);
            double dz =
                periodic_displacement(nj.pos.z() - p_i.pos.z(), domain_size);
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 < p_i.h * p_i.h && r2 > 1e-24) {
                double r = std::sqrt(r2);
                double dphi_dr_dummy, dW_dr;
                Kernels::adaptive_gravity_terms(r, p_i.h, dphi_dr_dummy, dW_dr);

                double V_j = nj.mass / nj.rho;
                Eigen::Vector3d dx_vec(dx, dy, dz);
                Eigen::Vector3d grad_W_i = -dW_dr * dx_vec / r;

                out.grad_rho += V_j * (nj.rho - p_i.rho) * grad_W_i;
                out.grad_p += V_j * (nj.pressure - p_i.pressure) * grad_W_i;
                out.grad_vx += V_j * (nj.vel.x() - p_i.vel.x()) * grad_W_i;
                out.grad_vy += V_j * (nj.vel.y() - p_i.vel.y()) * grad_W_i;
                out.grad_vz += V_j * (nj.vel.z() - p_i.vel.z()) * grad_W_i;
            }
        }
    }

    return out;
}

ReconstructedFace compute_face_reconstruction(
    const ParticleState& p_i, const ParticleGradients& grad_i,
    const ParticleState& p_j, const ParticleGradients& grad_j,
    double domain_size, double density_floor, double pressure_floor) {
    ReconstructedFace face;
    face.is_valid = false;

    double dx = periodic_displacement(p_j.pos.x() - p_i.pos.x(), domain_size);
    double dy = periodic_displacement(p_j.pos.y() - p_i.pos.y(), domain_size);
    double dz = periodic_displacement(p_j.pos.z() - p_i.pos.z(), domain_size);

    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < 1e-24) return face;

    face.r = std::sqrt(r2);
    Eigen::Vector3d dx_vec(dx, dy, dz);
    face.n = dx_vec / face.r;

    double fraction_i = p_i.h / (p_i.h + p_j.h);
    double fraction_j = 1.0 - fraction_i;

    Eigen::Vector3d dx_face_i = fraction_i * dx_vec;
    Eigen::Vector3d dx_face_j = -fraction_j * dx_vec;

#ifndef DISABLE_LIMITER
#ifdef THEORETICAL_LIMITER
    double rho_bar = p_i.rho + fraction_i * (p_j.rho - p_i.rho);
    double p_bar = p_i.pressure + fraction_i * (p_j.pressure - p_i.pressure);
    double vx_bar = p_i.vel.x() + fraction_i * (p_j.vel.x() - p_i.vel.x());
    double vy_bar = p_i.vel.y() + fraction_i * (p_j.vel.y() - p_i.vel.y());
    double vz_bar = p_i.vel.z() + fraction_i * (p_j.vel.z() - p_i.vel.z());

    face.rho_L = apply_pairwise_limiter(
        p_i.rho, p_j.rho, p_i.rho + grad_i.grad_rho.dot(dx_face_i), rho_bar);
    face.rho_R = apply_pairwise_limiter(
        p_j.rho, p_i.rho, p_j.rho + grad_j.grad_rho.dot(dx_face_j), rho_bar);

    face.p_L = apply_pairwise_limiter(
        p_i.pressure, p_j.pressure, p_i.pressure + grad_i.grad_p.dot(dx_face_i),
        p_bar);
    face.p_R = apply_pairwise_limiter(
        p_j.pressure, p_i.pressure, p_j.pressure + grad_j.grad_p.dot(dx_face_j),
        p_bar);

    face.v_L.x() = apply_pairwise_limiter(
        p_i.vel.x(), p_j.vel.x(), p_i.vel.x() + grad_i.grad_vx.dot(dx_face_i),
        vx_bar);
    face.v_R.x() = apply_pairwise_limiter(
        p_j.vel.x(), p_i.vel.x(), p_j.vel.x() + grad_j.grad_vx.dot(dx_face_j),
        vx_bar);

    face.v_L.y() = apply_pairwise_limiter(
        p_i.vel.y(), p_j.vel.y(), p_i.vel.y() + grad_i.grad_vy.dot(dx_face_i),
        vy_bar);
    face.v_R.y() = apply_pairwise_limiter(
        p_j.vel.y(), p_i.vel.y(), p_j.vel.y() + grad_j.grad_vy.dot(dx_face_j),
        vy_bar);

    face.v_L.z() = apply_pairwise_limiter(
        p_i.vel.z(), p_j.vel.z(), p_i.vel.z() + grad_i.grad_vz.dot(dx_face_i),
        vz_bar);
    face.v_R.z() = apply_pairwise_limiter(
        p_j.vel.z(), p_i.vel.z(), p_j.vel.z() + grad_j.grad_vz.dot(dx_face_j),
        vz_bar);

    face.rho_L = std::max(face.rho_L, density_floor);
    face.rho_R = std::max(face.rho_R, density_floor);
    face.p_L = std::max(face.p_L, pressure_floor);
    face.p_R = std::max(face.p_R, pressure_floor);

#else
    face.rho_L =
        std::max(p_i.rho + grad_i.grad_rho.dot(dx_face_i), density_floor);
    face.rho_R =
        std::max(p_j.rho + grad_j.grad_rho.dot(dx_face_j), density_floor);

    face.p_L =
        std::max(p_i.pressure + grad_i.grad_p.dot(dx_face_i), pressure_floor);
    face.p_R =
        std::max(p_j.pressure + grad_j.grad_p.dot(dx_face_j), pressure_floor);

    face.v_L.x() = p_i.vel.x() + grad_i.grad_vx.dot(dx_face_i);
    face.v_L.y() = p_i.vel.y() + grad_i.grad_vy.dot(dx_face_i);
    face.v_L.z() = p_i.vel.z() + grad_i.grad_vz.dot(dx_face_i);

    face.v_R.x() = p_j.vel.x() + grad_j.grad_vx.dot(dx_face_j);
    face.v_R.y() = p_j.vel.y() + grad_j.grad_vy.dot(dx_face_j);
    face.v_R.z() = p_j.vel.z() + grad_j.grad_vz.dot(dx_face_j);
#endif
#else
#ifdef ZEROTH_ORDER_RECONSTRUCTION
    face.rho_L = p_i.rho;
    face.rho_R = p_j.rho;

    face.p_L = p_i.pressure;
    face.p_R = p_j.pressure;

    face.v_L = p_i.vel;
    face.v_R = p_j.vel;

#else
    face.rho_L =
        std::max(p_i.rho + grad_i.grad_rho.dot(dx_face_i), 0.01 * p_i.rho);
    face.rho_R =
        std::max(p_j.rho + grad_j.grad_rho.dot(dx_face_j), 0.01 * p_j.rho);

    face.p_L = std::max(p_i.pressure + grad_i.grad_p.dot(dx_face_i),
                        0.01 * p_i.pressure);
    face.p_R = std::max(p_j.pressure + grad_j.grad_p.dot(dx_face_j),
                        0.01 * p_j.pressure);

    face.v_L.x() = p_i.vel.x() + grad_i.grad_vx.dot(dx_face_i);
    face.v_L.y() = p_i.vel.y() + grad_i.grad_vy.dot(dx_face_i);
    face.v_L.z() = p_i.vel.z() + grad_i.grad_vz.dot(dx_face_i);

    face.v_R.x() = p_j.vel.x() + grad_j.grad_vx.dot(dx_face_j);
    face.v_R.y() = p_j.vel.y() + grad_j.grad_vy.dot(dx_face_j);
    face.v_R.z() = p_j.vel.z() + grad_j.grad_vz.dot(dx_face_j);
#endif
#endif

    face.is_valid = true;
    return face;
}

}  // namespace Reconstruction