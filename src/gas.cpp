#include "gas.h"

#include <omp.h>

RiemannSolver::RiemannSolver(int mesh_size)
    : rho_L(mesh_size),
      p_L(mesh_size),
      vn_L(mesh_size),
      vt1_L(mesh_size),
      vt2_L(mesh_size),
      E_L(mesh_size),
      mom_n_L(mesh_size),
      mom_t1_L(mesh_size),
      mom_t2_L(mesh_size),
      rho_R(mesh_size),
      p_R(mesh_size),
      vn_R(mesh_size),
      vt1_R(mesh_size),
      vt2_R(mesh_size),
      E_R(mesh_size),
      mom_n_R(mesh_size),
      mom_t1_R(mesh_size),
      mom_t2_R(mesh_size),
      cs_L(mesh_size),
      cs_R(mesh_size),
      S_L(mesh_size),
      S_R(mesh_size),
      F_dens_L(mesh_size),
      F_dens_R(mesh_size),
      F_momn_L(mesh_size),
      F_momn_R(mesh_size),
      F_momt1_L(mesh_size),
      F_momt1_R(mesh_size),
      F_momt2_L(mesh_size),
      F_momt2_R(mesh_size),
      F_en_L(mesh_size),
      F_en_R(mesh_size),
      flux_density(mesh_size),
      flux_mom_n(mesh_size),
      flux_mom_t1(mesh_size),
      flux_mom_t2(mesh_size),
      flux_energy(mesh_size),
      ie_L(mesh_size),
      ie_R(mesh_size),
      F_ie_L(mesh_size),
      F_ie_R(mesh_size),
      flux_ie(mesh_size),
      flux_ie_sh(mesh_size),
      q_minus(mesh_size),
      q_plus(mesh_size),
      dq_L(mesh_size),
      dq_R(mesh_size),
      slope(mesh_size) {
    int total_cells = mesh_size * mesh_size * mesh_size;
    den_star.resize(total_cells);
    S_star.resize(total_cells);
    omega_L.resize(total_cells);
    omega_R.resize(total_cells);
    denom_L.resize(total_cells);
    denom_R.resize(total_cells);
    mom_n_star_L.resize(total_cells);
    mom_n_star_R.resize(total_cells);
    mom_t1_star_L.resize(total_cells);
    mom_t1_star_R.resize(total_cells);
    mom_t2_star_L.resize(total_cells);
    mom_t2_star_R.resize(total_cells);
    E_star_L.resize(total_cells);
    E_star_R.resize(total_cells);
    ie_star_L.resize(total_cells);
    ie_star_R.resize(total_cells);
}

inline void RiemannSolver::reconstruct_muscl(const Grid3D& q, Grid3D& q_L,
                                             Grid3D& q_R, int axis) {
    const int rollDir = -1;  // +1 cell shift (i+1)
    const int rollLeft = 1;  // -1 cell shift (i-1)

    q_minus = q.roll(rollLeft, axis);
    q_plus = q.roll(rollDir, axis);

    dq_L.data = q.array() - q_minus.array();
    dq_R.data = q_plus.array() - q.array();

    // Minmod logic done in-place to avoid returning new Grid3D
    slope.data = (dq_L.array() * dq_R.array() > 0.0)
                     .select(dq_L.array().abs().min(dq_R.array().abs()) *
                                 dq_L.array().sign(),
                             0.0);

    q_L.data = q.array() + 0.5 * slope.array();
    q_R.data = q_plus.array() - 0.5 * slope.roll(rollDir, axis).array();
}

inline void RiemannSolver::apply_hllc_flux(const Grid3D& F_L, const Grid3D& F_R,
                                           const Grid3D& U_L, const Grid3D& U_R,
                                           const Eigen::ArrayXd& U_star_L,
                                           const Eigen::ArrayXd& U_star_R,
                                           const Eigen::ArrayXd& S_star,
                                           Grid3D& F_out) {
    Eigen::ArrayXd F_star_L =
        F_L.array() + S_L.array() * (U_star_L - U_L.array());
    Eigen::ArrayXd F_star_R =
        F_R.array() + S_R.array() * (U_star_R - U_R.array());

    F_out.data =
        (S_L.array() >= 0.0)
            .select(
                F_L.array(),
                (S_star >= 0.0)
                    .select(
                        F_star_L,
                        (S_R.array() <= 0.0).select(F_R.array(), F_star_R)));
}

void RiemannSolver::compute_fluxes(const GasGrid& grid, int axis,
                                   double gamma) {
    const Grid3D& density = grid.get_density();
    const Grid3D& energy = grid.get_energy();
    const Grid3D& pressure = grid.get_pressure();
    const Grid3D& ie = grid.get_internal_energy();

    // References mapped to the Local Coordinate System based on the sweep axis
    const Grid3D& mom_n = (axis == 0) ? grid.get_momentum_x()
                                      : ((axis == 1) ? grid.get_momentum_y()
                                                     : grid.get_momentum_z());
    const Grid3D& v_n = (axis == 0) ? grid.get_velocity_x()
                                    : ((axis == 1) ? grid.get_velocity_y()
                                                   : grid.get_velocity_z());

    const Grid3D& mom_t1 = (axis == 0) ? grid.get_momentum_y()
                                       : ((axis == 1) ? grid.get_momentum_z()
                                                      : grid.get_momentum_x());
    const Grid3D& v_t1 = (axis == 0) ? grid.get_velocity_y()
                                     : ((axis == 1) ? grid.get_velocity_z()
                                                    : grid.get_velocity_x());

    const Grid3D& mom_t2 = (axis == 0) ? grid.get_momentum_z()
                                       : ((axis == 1) ? grid.get_momentum_x()
                                                      : grid.get_momentum_y());
    const Grid3D& v_t2 = (axis == 0) ? grid.get_velocity_z()
                                     : ((axis == 1) ? grid.get_velocity_x()
                                                    : grid.get_velocity_y());

    // Reconstruct primitive variables
    reconstruct_muscl(density, rho_L, rho_R, axis);
    reconstruct_muscl(pressure, p_L, p_R, axis);
    reconstruct_muscl(v_n, vn_L, vn_R, axis);
    reconstruct_muscl(v_t1, vt1_L, vt1_R, axis);
    reconstruct_muscl(v_t2, vt2_L, vt2_R, axis);
    reconstruct_muscl(ie, ie_L, ie_R, axis);

    // Reconstruct Energy and Momenta from the limited primitive variables
    // Doing this guarantees thermodynamic consistency at the interfaces
    E_L.data = p_L.array() / (gamma - 1.0) +
               0.5 * rho_L.array() *
                   (vn_L.array().square() + vt1_L.array().square() +
                    vt2_L.array().square());
    E_R.data = p_R.array() / (gamma - 1.0) +
               0.5 * rho_R.array() *
                   (vn_R.array().square() + vt1_R.array().square() +
                    vt2_R.array().square());

    mom_n_L.data = rho_L.array() * vn_L.array();
    mom_n_R.data = rho_R.array() * vn_R.array();
    mom_t1_L.data = rho_L.array() * vt1_L.array();
    mom_t1_R.data = rho_R.array() * vt1_R.array();
    mom_t2_L.data = rho_L.array() * vt2_L.array();
    mom_t2_R.data = rho_R.array() * vt2_R.array();

    // Compute Sound Speeds and Wave Signal Speeds
    cs_L.data = (gamma * p_L.array() / rho_L.array()).sqrt();
    cs_R.data = (gamma * p_R.array() / rho_R.array()).sqrt();

    // Stability: Zero out sound speed where density is effectively zero or NaN
    cs_L.data = (rho_L.array() > 1e-12).select(cs_L.array(), 0.0);
    cs_R.data = (rho_R.array() > 1e-12).select(cs_R.array(), 0.0);
    cs_L.data = (cs_L.array() == cs_L.array()).select(cs_L.array(), 0.0);
    cs_R.data = (cs_R.array() == cs_R.array()).select(cs_R.array(), 0.0);

    // HLL wave speed estimates (Davis 1988)
    S_L.data =
        (vn_L.array() - cs_L.array()).cwiseMin(vn_R.array() - cs_R.array());
    S_R.data =
        (vn_L.array() + cs_L.array()).cwiseMax(vn_R.array() + cs_R.array());

    // Evaluate Physical Fluxes for Left and Right states
    F_dens_L.data = rho_L.array() * vn_L.array();
    F_dens_R.data = rho_R.array() * vn_R.array();

    F_momn_L.data = rho_L.array() * vn_L.array().square() + p_L.array();
    F_momn_R.data = rho_R.array() * vn_R.array().square() + p_R.array();

    F_momt1_L.data = rho_L.array() * vn_L.array() * vt1_L.array();
    F_momt1_R.data = rho_R.array() * vn_R.array() * vt1_R.array();

    F_momt2_L.data = rho_L.array() * vn_L.array() * vt2_L.array();
    F_momt2_R.data = rho_R.array() * vn_R.array() * vt2_R.array();

    F_en_L.data = (E_L.array() + p_L.array()) * vn_L.array();
    F_en_R.data = (E_R.array() + p_R.array()) * vn_R.array();

    F_ie_L.data = ie_L.array() * vn_L.array();
    F_ie_R.data = ie_R.array() * vn_R.array();

    // --- THE HLLC RIEMANN SOLVER ---

    // Calculate the middle wave speed (S_star)
    den_star = rho_L.array() * (S_L.array() - vn_L.array()) -
               rho_R.array() * (S_R.array() - vn_R.array());
    den_star = (den_star.abs() < 1e-12).select(-1e-12, den_star);

    S_star = (p_R.array() - p_L.array() +
              rho_L.array() * vn_L.array() * (S_L.array() - vn_L.array()) -
              rho_R.array() * vn_R.array() * (S_R.array() - vn_R.array())) /
             den_star;

    // Calculate scaling factors for the "Star" regions
    omega_L = rho_L.array() * (S_L.array() - vn_L.array()) /
              (S_L.array() - S_star - 1e-12);
    omega_R = rho_R.array() * (S_R.array() - vn_R.array()) /
              (S_R.array() - S_star + 1e-12);

    // Prevent division by zero in the energy term
    denom_L = rho_L.array() * (S_L.array() - vn_L.array());
    denom_L = (denom_L.abs() < 1e-12).select(-1e-12, denom_L);
    denom_R = rho_R.array() * (S_R.array() - vn_R.array());
    denom_R = (denom_R.abs() < 1e-12).select(1e-12, denom_R);

    // Compute the Star States (U_star)
    const Eigen::ArrayXd& rho_star_L = omega_L;
    const Eigen::ArrayXd& rho_star_R = omega_R;

    mom_n_star_L = omega_L * S_star;
    mom_n_star_R = omega_R * S_star;

    mom_t1_star_L = omega_L * vt1_L.array();
    mom_t1_star_R = omega_R * vt1_R.array();

    mom_t2_star_L = omega_L * vt2_L.array();
    mom_t2_star_R = omega_R * vt2_R.array();

    E_star_L =
        omega_L * (E_L.array() / rho_L.array() +
                   (S_star - vn_L.array()) * (S_star + p_L.array() / denom_L));
    E_star_R =
        omega_R * (E_R.array() / rho_R.array() +
                   (S_star - vn_R.array()) * (S_star + p_R.array() / denom_R));

    ie_star_L = omega_L * (ie_L.array() / rho_L.array());
    ie_star_R = omega_R * (ie_R.array() / rho_R.array());

    // Calculate all intercell fluxes
    apply_hllc_flux(F_dens_L, F_dens_R, rho_L, rho_R, rho_star_L, rho_star_R,
                    S_star, flux_density);
    apply_hllc_flux(F_momn_L, F_momn_R, mom_n_L, mom_n_R, mom_n_star_L,
                    mom_n_star_R, S_star, flux_mom_n);
    apply_hllc_flux(F_momt1_L, F_momt1_R, mom_t1_L, mom_t1_R, mom_t1_star_L,
                    mom_t1_star_R, S_star, flux_mom_t1);
    apply_hllc_flux(F_momt2_L, F_momt2_R, mom_t2_L, mom_t2_R, mom_t2_star_L,
                    mom_t2_star_R, S_star, flux_mom_t2);
    apply_hllc_flux(F_en_L, F_en_R, E_L, E_R, E_star_L, E_star_R, S_star,
                    flux_energy);
    apply_hllc_flux(F_ie_L, F_ie_R, ie_L, ie_R, ie_star_L, ie_star_R, S_star,
                    flux_ie);

    // Pre-calculate shifted fluxes for the grid update
    flux_density_sh = flux_density.roll(1, axis);
    flux_mom_n_sh = flux_mom_n.roll(1, axis);
    flux_mom_t1_sh = flux_mom_t1.roll(1, axis);
    flux_mom_t2_sh = flux_mom_t2.roll(1, axis);
    flux_energy_sh = flux_energy.roll(1, axis);
    flux_ie_sh = flux_ie.roll(1, axis);
}

GasGrid::GasGrid(const Config& conf)
    : density(conf.MESH_SIZE),
      momentum_x(conf.MESH_SIZE),
      momentum_y(conf.MESH_SIZE),
      momentum_z(conf.MESH_SIZE),
      energy(conf.MESH_SIZE),
      pressure(conf.MESH_SIZE),
      velocity_x(conf.MESH_SIZE),
      velocity_y(conf.MESH_SIZE),
      velocity_z(conf.MESH_SIZE),
      internal_energy(conf.MESH_SIZE),
      solver(conf.MESH_SIZE),
      config(conf) {}

void GasGrid::update_primitive_variables() {
    const int total_cells =
        config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;

    // Pre-calculate constants to avoid doing it inside the loop
    const double gamma_minus_1 = config.GAMMA - 1.0;
    const double inv_gamma_minus_1 = 1.0 / gamma_minus_1;
    const double floor = 1e-12;

    // Extract raw pointers so the compiler can use vectorization
    double* rho = density.array().data();
    double* mx = momentum_x.array().data();
    double* my = momentum_y.array().data();
    double* mz = momentum_z.array().data();
    double* en = energy.array().data();
    double* vx = velocity_x.array().data();
    double* vy = velocity_y.array().data();
    double* vz = velocity_z.array().data();
    double* pr = pressure.array().data();
    double* ie_tracked = internal_energy.array().data();

#pragma omp parallel for schedule(static)
    for (int i = 0; i < total_cells; ++i) {
        double local_rho = rho[i];

        if (local_rho > floor) {
            double inv_rho = 1.0 / local_rho;
            double local_vx = mx[i] * inv_rho;
            double local_vy = my[i] * inv_rho;
            double local_vz = mz[i] * inv_rho;

            vx[i] = local_vx;
            vy[i] = local_vy;
            vz[i] = local_vz;

            double ke = 0.5 * local_rho *
                        (local_vx * local_vx + local_vy * local_vy +
                         local_vz * local_vz);
            double e_tot = en[i];
            double p;

            // --- THE DUAL ENERGY SWITCH ---
            // If internal energy is > 0.1% of total energy, trust the total
            // energy
            if (e_tot > 0.0 && (e_tot - ke) / e_tot > 0.001) {
                p = gamma_minus_1 * (e_tot - ke);
                ie_tracked[i] =
                    p * inv_gamma_minus_1;  // Sync tracked IE to prevent drift
            } else {
                // Hypersonic regime! Total energy is polluted. Trust tracked
                // IE.
                p = gamma_minus_1 * ie_tracked[i];
                en[i] = p * inv_gamma_minus_1 + ke;  // Sync Total Energy
            }

            if (p < floor) {
                p = floor;
                ie_tracked[i] = floor * inv_gamma_minus_1;
                en[i] = ie_tracked[i] + ke;
            }

            pr[i] = p;
        } else {
            // Vacuum cell handling
            vx[i] = 0.0;
            vy[i] = 0.0;
            vz[i] = 0.0;
            pr[i] = floor;
            ie_tracked[i] = floor * inv_gamma_minus_1;
            en[i] = ie_tracked[i];
        }
    }
}

void GasGrid::compute_and_apply_fluxes(double dt) {
    double factor = dt / config.CELL_SIZE;

    for (int axis = 0; axis < 3; ++axis) {
        update_primitive_variables();
        solver.compute_fluxes(*this, axis, config.GAMMA);

        // Helper references to update original coordinates
        Grid3D* m_n;
        Grid3D* m_t1;
        Grid3D* m_t2;
        Grid3D* v_n;
        if (axis == 0) {
            m_n = &momentum_x;
            m_t1 = &momentum_y;
            m_t2 = &momentum_z;
            v_n = &velocity_x;
        } else if (axis == 1) {
            m_n = &momentum_y;
            m_t1 = &momentum_z;
            m_t2 = &momentum_x;
            v_n = &velocity_y;
        } else {
            m_n = &momentum_z;
            m_t1 = &momentum_x;
            m_t2 = &momentum_y;
            v_n = &velocity_z;
        }

        // Central difference derivative for divergence: (v[i+1] - v[i-1]) / 2dx
        Grid3D dv(config.MESH_SIZE);
        dv.data = v_n->roll(-1, axis).array() - v_n->roll(1, axis).array();

        int total_cells =
            config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;

        double* d_rho = density.array().data();
        double* d_mn = m_n->array().data();
        double* d_mt1 = m_t1->array().data();
        double* d_mt2 = m_t2->array().data();
        double* d_en = energy.array().data();
        double* d_ie = internal_energy.array().data();

        double* p_arr = pressure.array().data();
        double* dv_arr = dv.array().data();

        const double* f_rho = solver.get_flux_density().array().data();
        const double* f_mn = solver.get_flux_mom_n().array().data();
        const double* f_mt1 = solver.get_flux_mom_t1().array().data();
        const double* f_mt2 = solver.get_flux_mom_t2().array().data();
        const double* f_en = solver.get_flux_energy().array().data();
        const double* f_ie = solver.get_flux_ie().array().data();

        const double* f_rho_sh = solver.get_flux_density_sh().array().data();
        const double* f_mn_sh = solver.get_flux_mom_n_sh().array().data();
        const double* f_mt1_sh = solver.get_flux_mom_t1_sh().array().data();
        const double* f_mt2_sh = solver.get_flux_mom_t2_sh().array().data();
        const double* f_en_sh = solver.get_flux_energy_sh().array().data();
        const double* f_ie_sh = solver.get_flux_ie_sh().array().data();

#pragma omp parallel for simd schedule(static)
        for (int i = 0; i < total_cells; ++i) {
            d_rho[i] -= factor * (f_rho[i] - f_rho_sh[i]);
            if (d_rho[i] < 1e-12) d_rho[i] = 1e-12;  // Density floor

            d_mn[i] -= factor * (f_mn[i] - f_mn_sh[i]);
            d_mt1[i] -= factor * (f_mt1[i] - f_mt1_sh[i]);
            d_mt2[i] -= factor * (f_mt2[i] - f_mt2_sh[i]);
            d_en[i] -= factor * (f_en[i] - f_en_sh[i]);

            d_ie[i] -= factor * (f_ie[i] - f_ie_sh[i]);
            d_ie[i] -= 0.5 * factor * p_arr[i] * dv_arr[i];  // PdV work
        }
    }
}

void GasGrid::hydro_step(double dt) {
    // Backup the current conservative state (U^n)
    Grid3D old_rho = density;
    Grid3D old_mx = momentum_x;
    Grid3D old_my = momentum_y;
    Grid3D old_mz = momentum_z;
    Grid3D old_en = energy;
    Grid3D old_ie = internal_energy;

    // RK2 Stage 1
    // Advance using the initial state: U^{(1)} = U^n + dt * L(U^n)
    compute_and_apply_fluxes(dt);

    // RK2 Stage 2
    // Advance again using the intermediate state: U^{(2)} = U^{(1)} + dt *
    // L(U^{(1)})
    compute_and_apply_fluxes(dt);

    // Final Averaging (FUSED LOOP)
    int total_cells = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;

    double* d_rho = density.array().data();
    double* d_mx = momentum_x.array().data();
    double* d_my = momentum_y.array().data();
    double* d_mz = momentum_z.array().data();
    double* d_en = energy.array().data();
    double* d_ie = internal_energy.array().data();

    double* o_rho = old_rho.array().data();
    double* o_mx = old_mx.array().data();
    double* o_my = old_my.array().data();
    double* o_mz = old_mz.array().data();
    double* o_en = old_en.array().data();
    double* o_ie = old_ie.array().data();

#pragma omp parallel for simd schedule(static)
    for (int i = 0; i < total_cells; ++i) {
        d_rho[i] = 0.5 * (o_rho[i] + d_rho[i]);
        d_mx[i] = 0.5 * (o_mx[i] + d_mx[i]);
        d_my[i] = 0.5 * (o_my[i] + d_my[i]);
        d_mz[i] = 0.5 * (o_mz[i] + d_mz[i]);
        d_en[i] = 0.5 * (o_en[i] + d_en[i]);
        d_ie[i] = 0.5 * (o_ie[i] + d_ie[i]);
    }

    // Update primitive variables once at the end so the gravity step has the
    // right values
    update_primitive_variables();
}

double GasGrid::get_cfl_timestep() const {
    if (!config.USE_HYDRO) return std::numeric_limits<double>::infinity();

    Grid3D cs_sq(config.MESH_SIZE);
    cs_sq.data = (config.GAMMA * pressure.array() / density.array());
    cs_sq.data = (density.array() > 1e-12).select(cs_sq.array(), 0.0);
    cs_sq.data = (cs_sq.array() == cs_sq.array()).select(cs_sq.array(), 0.0);

    Grid3D v_mag(config.MESH_SIZE);
    v_mag.data = (velocity_x.array().square() + velocity_y.array().square() +
                  velocity_z.array().square())
                     .sqrt();

    Grid3D signal_vel(config.MESH_SIZE);
    signal_vel.data = v_mag.array() + cs_sq.array().sqrt();

    double max_signal_vel = signal_vel.maxCoeff();
    if (max_signal_vel < 1e-9) return std::numeric_limits<double>::infinity();

    return (config.CELL_SIZE / max_signal_vel) * config.CFL_SAFETY_FACTOR;
}