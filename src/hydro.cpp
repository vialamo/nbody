#include "hydro.h"

#include <omp.h>

#include "config.h"
#include "debug_utils.h"
#include "kernels.h"
#include "math_utils.h"
#include "mfm.h"
#include "riemann.h"

namespace HydroSolver {

namespace {  // Anonymous namespace for internal hydro helpers

// Computes the maximum signal velocity between two particles and updates
// their cached maximums using lock-free atomics for the CFL timestep
// condition
static inline void update_signal_velocities(
    const Reconstruction::ParticleState& p_i,
    const Reconstruction::ParticleState& p_j, double dx, double dy, double dz,
    double r2, double a, double a_inv, double a_inv3, double gamma,
    double& v_sig_max_i, double& v_sig_max_j, double& out_v_sig_ij,
    double& out_c_phys_j) {
    double r = std::sqrt(r2);
    double rho_phys_i = p_i.rho * a_inv3;
    double p_phys_i = p_i.pressure * a_inv;
    double c_phys_i =
        (rho_phys_i > 1e-12) ? std::sqrt(gamma * p_phys_i / rho_phys_i) : 0.0;

    double rho_phys_j = p_j.rho * a_inv3;
    double p_phys_j = p_j.pressure * a_inv;
    double c_phys_j =
        (rho_phys_j > 1e-12) ? std::sqrt(gamma * p_phys_j / rho_phys_j) : 0.0;

    out_c_phys_j = c_phys_j;

    double dvx_phys = (p_j.vel.x() - p_i.vel.x()) * a;
    double dvy_phys = (p_j.vel.y() - p_i.vel.y()) * a;
    double dvz_phys = (p_j.vel.z() - p_i.vel.z()) * a;

    double dv_dot_dx_phys = dvx_phys * dx + dvy_phys * dy + dvz_phys * dz;
    double v_sig_ij_phys = c_phys_i + c_phys_j;

    if (dv_dot_dx_phys < 0.0) {
        v_sig_ij_phys -= dv_dot_dx_phys / r;
    }

    out_v_sig_ij = v_sig_ij_phys;

    // Fast atomic max using Compare-and-Swap (CAS) to avoid OpenMP lock
    // overhead
    auto atomic_update_max = [](double* address, double val) {
        uint64_t* addr_as_int = reinterpret_cast<uint64_t*>(address);
        uint64_t old_val = *addr_as_int;
        while (val > *reinterpret_cast<double*>(&old_val)) {
            uint64_t new_val_int = *reinterpret_cast<uint64_t*>(&val);
            if (__sync_bool_compare_and_swap(addr_as_int, old_val,
                                             new_val_int)) {
                break;
            }
            old_val = *addr_as_int;
        }
    };

    // Update particle i
    // Because of the (j <= i) symmetry optimization, thread 'j' might skip
    // 'i'. Therefore, thread 'i' MUST push the max velocity to 'j'
    atomic_update_max(&v_sig_max_i, v_sig_ij_phys);
    atomic_update_max(&v_sig_max_j, v_sig_ij_phys);
}

}  // end anonymous namespace

void compute_gradients(GasParticleSystem& gas, const Config& config) {
    if (gas.num_particles == 0) return;
    const Vec3& domain_size = config.domain_size;
    size_t step_ill_conditioned = 0;

    if (gas.num_active == 0) return;

    Reconstruction::FluidStateArrays system_data{
        gas.pos_x.data(),      gas.pos_y.data(),      gas.pos_z.data(),
        gas.vel_x_pred.data(), gas.vel_y_pred.data(), gas.vel_z_pred.data(),
        gas.mass.data(),       gas.rho.data(),        gas.pressure_pred.data(),
        gas.h.data(),          gas.cond_num.data()};

    bool debug_printed = false;

#pragma omp parallel reduction(+ : step_ill_conditioned)
    {
        // Declare the reusable buffer once per thread
        std::vector<int> neighbors;
        neighbors.reserve(256);

#pragma omp for schedule(dynamic, 64)
        for (size_t k = 0; k < gas.num_active; ++k) {
            size_t i = gas.active_indices[k];
            neighbors.clear();
            double local_max_rel_ke = 0.0;

            int stack[128];
            int stack_ptr = 0;
            stack[stack_ptr++] = 0;

            while (stack_ptr > 0) {
                int node_idx = stack[--stack_ptr];
                const BVHNode& node = gas.bvh_nodes[node_idx];

                double dist_sq =
                    min_periodic_dist_sq(gas.pos_x[i], node.bbox.min_x,
                                         node.bbox.max_x, domain_size.x) +
                    min_periodic_dist_sq(gas.pos_y[i], node.bbox.min_y,
                                         node.bbox.max_y, domain_size.y) +
                    min_periodic_dist_sq(gas.pos_z[i], node.bbox.min_z,
                                         node.bbox.max_z, domain_size.z);

                if (dist_sq > gas.h[i] * gas.h[i]) continue;

                if (node.particle_idx != -1) {
                    int j = node.particle_idx;

                    double dx = periodic_displacement(
                        gas.pos_x[j] - gas.pos_x[i], domain_size.x);
                    double dy = periodic_displacement(
                        gas.pos_y[j] - gas.pos_y[i], domain_size.y);
                    double dz = periodic_displacement(
                        gas.pos_z[j] - gas.pos_z[i], domain_size.z);
                    double r2 = dx * dx + dy * dy + dz * dz;

                    if ((r2 < gas.h[i] * gas.h[i] ||
                         r2 < gas.h[j] * gas.h[j]) &&
                        r2 > 1e-24) {
                        double rel_vx = gas.vel_x[j] - gas.vel_x[i];
                        double rel_vy = gas.vel_y[j] - gas.vel_y[i];
                        double rel_vz = gas.vel_z[j] - gas.vel_z[i];
                        double rel_v2 =
                            rel_vx * rel_vx + rel_vy * rel_vy + rel_vz * rel_vz;
                        local_max_rel_ke =
                            std::max(local_max_rel_ke, 0.5 * rel_v2);
                        if (r2 < gas.h[i] * gas.h[i]) {
                            neighbors.push_back(j);
                        }
                    }
                } else {
                    stack[stack_ptr++] = node.left_child;
                    stack[stack_ptr++] = node.right_child;
                }
            }

            gas.max_rel_ke[i] = local_max_rel_ke;

            // Calculate delta E_grav = |a_grav| * h
            double a_grav_mag = std::sqrt(gas.acc_x[i] * gas.acc_x[i] +
                                          gas.acc_y[i] * gas.acc_y[i] +
                                          gas.acc_z[i] * gas.acc_z[i]);
            gas.delta_E_grav[i] = a_grav_mag * gas.h[i];

            bool ill_conditioned =
                (gas.cond_num[i] < 0.0 || gas.cond_num[i] > N_cond_crit);

            Reconstruction::ParticleGradients grads =
                Reconstruction::compute_single_particle_gradients(
                    i, neighbors.data(), neighbors.size(), system_data,
                    gas.B_matrix[i], ill_conditioned, domain_size);

            // Map results back to SOA
            gas.grad_rho[i] = grads.grad_rho;
            gas.grad_p[i] = grads.grad_p;
            gas.grad_vx[i] = grads.grad_vx;
            gas.grad_vy[i] = grads.grad_vy;
            gas.grad_vz[i] = grads.grad_vz;

            if (ill_conditioned) {
                step_ill_conditioned++;
                // DebugUtils::print_ill_conditioned_diagnostic(debug_printed,
                // i,
                //    gas, neighbors, domain_size);
            }
            // DebugUtils::print_gradients(gas, i);
        }
    }
    gas.ill_conditioned_cases += step_ill_conditioned;
}

void compute_forces(GasParticleSystem& gas, const Config& config, double a) {
    if (gas.num_particles == 0) return;
    const Vec3& domain_size = config.domain_size;

    if (gas.num_active == 0) return;

    // Reset accumulators
    for (size_t k = 0; k < gas.num_active; ++k) {
        size_t i = gas.active_indices[k];
        gas.hydro_acc_x[i] = 0.0;
        gas.hydro_acc_y[i] = 0.0;
        gas.hydro_acc_z[i] = 0.0;
        gas.du_dt[i] = 0.0;
        gas.de_dt[i] = 0.0;
        gas.v_sig_max[i] = 0.0;
    }

    size_t step_sph_fallbacks = 0;
    size_t riemann_failed_cases = 0;

#pragma omp parallel for schedule(dynamic, 64) \
    reduction(+ : step_sph_fallbacks, riemann_failed_cases)
    for (size_t k = 0; k < gas.num_active; ++k) {
        size_t i = gas.active_indices[k];

        Reconstruction::ParticleState p_i;
        p_i.pos = Eigen::Vector3d(gas.pos_x[i], gas.pos_y[i], gas.pos_z[i]);
        p_i.mass = gas.mass[i];
        p_i.rho = gas.rho[i];
        p_i.h = gas.h[i];
        p_i.vel = Eigen::Vector3d(gas.vel_x_pred[i], gas.vel_y_pred[i],
                                  gas.vel_z_pred[i]);
        p_i.pressure = gas.pressure_pred[i];

        // Build the isolated state proxy for particle i
        Reconstruction::ParticleGradients grad_i;
        grad_i.grad_rho = gas.grad_rho[i];
        grad_i.grad_p = gas.grad_p[i];
        grad_i.grad_vx = gas.grad_vx[i];
        grad_i.grad_vy = gas.grad_vy[i];
        grad_i.grad_vz = gas.grad_vz[i];
        grad_i.ill_conditioned =
            (gas.cond_num[i] < 0.0 || gas.cond_num[i] > N_cond_crit);

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = gas.bvh_nodes[node_idx];

            double dist_sq =
                min_periodic_dist_sq(p_i.pos.x(), node.bbox.min_x,
                                     node.bbox.max_x, domain_size.x) +
                min_periodic_dist_sq(p_i.pos.y(), node.bbox.min_y,
                                     node.bbox.max_y, domain_size.y) +
                min_periodic_dist_sq(p_i.pos.z(), node.bbox.min_z,
                                     node.bbox.max_z, domain_size.z);

            if (dist_sq > p_i.h * p_i.h) continue;

            if (node.particle_idx != -1) {
                int j = node.particle_idx;

                // Maintain N3L conservation for active pairs, but ensure
                // inactive neighbors are always processed by the active
                // particle
                if (gas.is_active[j] && j <= static_cast<int>(i)) continue;

                double dx = periodic_displacement(gas.pos_x[j] - p_i.pos.x(),
                                                  domain_size.x);
                double dy = periodic_displacement(gas.pos_y[j] - p_i.pos.y(),
                                                  domain_size.y);
                double dz = periodic_displacement(gas.pos_z[j] - p_i.pos.z(),
                                                  domain_size.z);
                double r2 = dx * dx + dy * dy + dz * dz;

                if (r2 < 1e-24) continue;

                // Particles interact if they fall within either smoothing
                // length
                if (r2 < p_i.h * p_i.h || r2 < gas.h[j] * gas.h[j]) {
                    Reconstruction::ParticleState p_j;
                    p_j.pos = Eigen::Vector3d(gas.pos_x[j], gas.pos_y[j],
                                              gas.pos_z[j]);
                    p_j.mass = gas.mass[j];
                    p_j.rho = gas.rho[j];
                    p_j.h = gas.h[j];
                    p_j.vel =
                        Eigen::Vector3d(gas.vel_x_pred[j], gas.vel_y_pred[j],
                                        gas.vel_z_pred[j]);
                    p_j.pressure = gas.pressure_pred[j];

                    Reconstruction::ParticleGradients grad_j;
                    grad_j.grad_rho = gas.grad_rho[j];
                    grad_j.grad_p = gas.grad_p[j];
                    grad_j.grad_vx = gas.grad_vx[j];
                    grad_j.grad_vy = gas.grad_vy[j];
                    grad_j.grad_vz = gas.grad_vz[j];
                    grad_j.ill_conditioned = (gas.cond_num[j] < 0.0 ||
                                              gas.cond_num[j] > N_cond_crit);

                    // Solver pipeline
                    Reconstruction::ReconstructedFace face =
                        compute_face_reconstruction(
                            p_i, grad_i, p_j, grad_j, gas.B_matrix[i],
                            gas.B_matrix[j], domain_size, density_floor,
                            gas.pressure_floor);

                    if (!face.is_valid) continue;
                    if (face.used_sph_fallback) step_sph_fallbacks++;

#ifdef USE_MIDPOINT_QUADRATURE
                    double fraction_i = 0.5;
#else
                    double fraction_i = p_i.h / (p_i.h + p_j.h);
#endif
                    Eigen::Vector3d v_frame =
                        p_i.vel + fraction_i * (p_j.vel - p_i.vel);

                    // Convert to physical units
                    double a_inv = 1.0 / a;
                    double a_inv3 = a_inv * a_inv * a_inv;

                    Reconstruction::ReconstructedFace face_phys = face;
                    face_phys.rho_L = face.rho_L * a_inv3;
                    face_phys.rho_R = face.rho_R * a_inv3;
                    face_phys.p_L = face.p_L * a_inv;
                    face_phys.p_R = face.p_R * a_inv;
                    face_phys.v_L = face.v_L * a;
                    face_phys.v_R = face.v_R * a;

                    Eigen::Vector3d v_frame_phys = v_frame * a;

                    Riemann::FaceFlux flux_1d_phys;
                    if (config.use_exact_riemann) {
                        flux_1d_phys = Riemann::solve_exact(
                            face_phys, v_frame_phys, config.gamma);
                    } else {
                        flux_1d_phys = Riemann::solve_hllc(
                            face_phys, v_frame_phys, config.gamma);
                    }

                    // Fallback to First-Order Reconstruction if solvers still
                    // fail
                    if (!flux_1d_phys.success) {
                        // Overwrite the extrapolated face states with the
                        // particle center values
                        face_phys.rho_L = p_i.rho * a_inv3;
                        face_phys.rho_R = p_j.rho * a_inv3;
                        face_phys.p_L = p_i.pressure * a_inv;
                        face_phys.p_R = p_j.pressure * a_inv;
                        face_phys.v_L = p_i.vel * a;
                        face_phys.v_R = p_j.vel * a;

                        // Retry the exact solver with the flat states
                        flux_1d_phys = Riemann::solve_exact(
                            face_phys, v_frame_phys, config.gamma);
                    }

                    // Final fallback: Zero-Flux State
                    if (!flux_1d_phys.success) {
                        DebugUtils::print_riemann_diagnostic(
                            i, j, face_phys, v_frame_phys, config.gamma,
                            flux_1d_phys.error);
                        flux_1d_phys.P_star = 0.0;
                        flux_1d_phys.S_star = 0.0;
                        flux_1d_phys.flux_mom = Eigen::Vector3d::Zero();
                        riemann_failed_cases++;
                    }

                    // Calculate Maximum Signal Velocity for the CFL
                    // condition
                    double v_sig_ij, c_phys_j;
                    update_signal_velocities(p_i, p_j, dx, dy, dz, r2, a, a_inv,
                                             a_inv3, config.gamma,
                                             gas.v_sig_max[i], gas.v_sig_max[j],
                                             v_sig_ij, c_phys_j);

                    // Back to code units
                    double P_star_com = flux_1d_phys.P_star * a;
                    double S_star_com = flux_1d_phys.S_star * a_inv;

                    Eigen::Vector3d Force_mom = P_star_com * face.area_vec;
                    Eigen::Vector3d v_star_lab =
                        v_frame + (S_star_com * face.n);

                    double work_i =
                        P_star_com * (v_star_lab - p_i.vel).dot(face.area_vec);
                    double du_dt_i = -work_i / p_i.mass;

                    double Rate_energy =
                        P_star_com * v_star_lab.dot(face.area_vec);
                    double de_dt_i = -Rate_energy / p_i.mass;

                    // Apply the flux to the active particle 'i'
#pragma omp atomic
                    gas.hydro_acc_x[i] -= Force_mom.x() / p_i.mass;
#pragma omp atomic
                    gas.hydro_acc_y[i] -= Force_mom.y() / p_i.mass;
#pragma omp atomic
                    gas.hydro_acc_z[i] -= Force_mom.z() / p_i.mass;
#pragma omp atomic
                    gas.du_dt[i] += du_dt_i;
#pragma omp atomic
                    gas.de_dt[i] += de_dt_i;

                    // Apply the reverse flux to j if j is ALSO active
                    if (gas.is_active[j]) {
                        double work_j =
                            P_star_com *
                            (v_star_lab - p_j.vel).dot(face.area_vec);
                        double du_dt_j = work_j / p_j.mass;
                        double de_dt_j = Rate_energy / p_j.mass;
#pragma omp atomic
                        gas.hydro_acc_x[j] += Force_mom.x() / p_j.mass;
#pragma omp atomic
                        gas.hydro_acc_y[j] += Force_mom.y() / p_j.mass;
#pragma omp atomic
                        gas.hydro_acc_z[j] += Force_mom.z() / p_j.mass;
#pragma omp atomic
                        gas.du_dt[j] += du_dt_j;
#pragma omp atomic
                        gas.de_dt[j] += de_dt_j;
                    } else if (config.individual_particle_timesteps &&
                               !gas.needs_wakeup[j]) {
                        // OpenGadget3 Wake-up Limiter
                        // Wake up the sleeping neighbor if the signal
                        // velocity of the interaction exceeds 3x its own
                        // sound speed (f_w = 3.0)
                        bool og3_trigger = (v_sig_ij > 3.0 * c_phys_j);

                        // Saitoh & Makino (2009) Timestep Limiter
                        // Wake up the sleeping neighbor if its timestep is
                        // larger (e.g., 4x) than the active particle
                        // hitting it
                        // bool gizmo_trigger = (dt_step[j] > 4.0 *
                        // dt_step[i]);
                        if (og3_trigger /*|| gizmo_trigger*/) {
#pragma omp atomic write
                            gas.needs_wakeup[j] = 1;
                        }
                    }
                }
            } else {
                stack[stack_ptr++] = node.left_child;
                stack[stack_ptr++] = node.right_child;
            }
        }
    }
    gas.sph_fallback_cases += step_sph_fallbacks;
    gas.riemann_failed_cases += riemann_failed_cases;

    // For sod test, keep boundaries static
    if (config.initial_setup == InitialSetup::SodShockTube) {
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < gas.num_particles; ++i) {
            if (gas.pos_x[i] < 0.05 * config.domain_size.x ||
                gas.pos_x[i] > 0.95 * config.domain_size.x) {
                gas.hydro_acc_x[i] = 0.0;
                gas.hydro_acc_y[i] = 0.0;
                gas.hydro_acc_z[i] = 0.0;
                gas.du_dt[i] = 0.0;
            }
        }
    }
}

}  // namespace HydroSolver
