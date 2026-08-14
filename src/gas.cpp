#include "gas.h"

#include <omp.h>

#include <iostream>

#include "cmath"
#include "cooling.h"

constexpr double density_floor = 1e-12;

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
      Z(mesh_size),
      Z_L(mesh_size),
      Z_R(mesh_size),
      ie_L(mesh_size),
      ie_R(mesh_size),
      F_ie_L(mesh_size),
      F_ie_R(mesh_size),
      flux_ie(mesh_size),
      flux_ie_sh(mesh_size),
      flux_metal(mesh_size),
      flux_metal_sh(mesh_size),
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
    constexpr int rollRight = -1;  // +1 cell shift (i+1)
    constexpr int rollLeft = 1;    // -1 cell shift (i-1)

    q_minus = q.roll(rollLeft, axis);
    q_plus = q.roll(rollRight, axis);

    dq_L.data = q.array() - q_minus.array();
    dq_R.data = q_plus.array() - q.array();

    // Minmod logic done in-place
    slope.data = (dq_L.array() * dq_R.array() > 0.0)
                     .select(dq_L.array().abs().min(dq_R.array().abs()) *
                                 dq_L.array().sign(),
                             0.0);

    q_L.data = q.array() + 0.5 * slope.array();
    q_R.data = q_plus.array() - 0.5 * slope.roll(rollRight, axis).array();
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
    const Grid3D& pressure = grid.get_pressure();
    const Grid3D& ie = grid.get_internal_energy();
    const Grid3D& metal_dens = grid.get_metal_density();

    // Important: Advect conserved variables,
    // but reconstruct primitive variables
    Z.data = metal_dens.array() / density.array();

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
    reconstruct_muscl(Z, Z_L, Z_R, axis);

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
    cs_L.data = (rho_L.array() > density_floor).select(cs_L.array(), 0.0);
    cs_R.data = (rho_R.array() > density_floor).select(cs_R.array(), 0.0);
    cs_L.data = (cs_L.array() == cs_L.array()).select(cs_L.array(), 0.0);
    cs_R.data = (cs_R.array() == cs_R.array()).select(cs_R.array(), 0.0);

    // HLL wave speed estimates
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

    // HLLC RIEMANN SOLVER

    // Calculate the middle wave speed (S_star)
    constexpr double epsilon = 1e-15;
    den_star = rho_L.array() * (S_L.array() - vn_L.array()) -
               rho_R.array() * (S_R.array() - vn_R.array());

    den_star = (den_star.abs() < epsilon).select(-epsilon, den_star);

    S_star = (p_R.array() - p_L.array() +
              rho_L.array() * vn_L.array() * (S_L.array() - vn_L.array()) -
              rho_R.array() * vn_R.array() * (S_R.array() - vn_R.array())) /
             den_star;

    // Calculate scaling factors for the "Star" regions
    omega_L = rho_L.array() * (S_L.array() - vn_L.array()) /
              (S_L.array() - S_star - density_floor);
    omega_R = rho_R.array() * (S_R.array() - vn_R.array()) /
              (S_R.array() - S_star + density_floor);

    // Prevent division by zero in the energy term
    denom_L = rho_L.array() * (S_L.array() - vn_L.array());
    denom_L = (denom_L.abs() < epsilon).select(-epsilon, denom_L);
    denom_R = rho_R.array() * (S_R.array() - vn_R.array());
    denom_R = (denom_R.abs() < epsilon).select(epsilon, denom_R);

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
    flux_metal.data = (S_star >= 0.0)
                          .select(flux_density.array() * Z_L.array(),
                                  flux_density.array() * Z_R.array());

    // Pre-calculate shifted fluxes for the grid update
    flux_density_sh = flux_density.roll(1, axis);
    flux_mom_n_sh = flux_mom_n.roll(1, axis);
    flux_mom_t1_sh = flux_mom_t1.roll(1, axis);
    flux_mom_t2_sh = flux_mom_t2.roll(1, axis);
    flux_energy_sh = flux_energy.roll(1, axis);
    flux_ie_sh = flux_ie.roll(1, axis);
    flux_metal_sh = flux_metal.roll(1, axis);
}

GasGrid::GasGrid(const Config& conf)
    : density(conf.mesh_size),
      momentum_x(conf.mesh_size),
      momentum_y(conf.mesh_size),
      momentum_z(conf.mesh_size),
      energy(conf.mesh_size),
      pressure(conf.mesh_size),
      velocity_x(conf.mesh_size),
      velocity_y(conf.mesh_size),
      velocity_z(conf.mesh_size),
      metal_density(conf.mesh_size),
      internal_energy(conf.mesh_size),
      solver(conf.mesh_size),
      cooling_failed_cells(0),
      cooling_total_cycles(0),
      accumulated_radiated_energy(0.0),
      accumulated_photoheating_energy(0.0),
      accumulated_gravitational_work(0.0),
      accumulated_expansion_work(0.0),
      pressure_floor(0.0),
      config(conf) {
    double u_code_1K = Cooling::get_internal_energy_from_temp(1.0, 1.0, config);
    pressure_floor = (config.gamma - 1.0) * density_floor * u_code_1K;
}

void GasGrid::update_primitive_variables() {
    const int total_cells =
        config.mesh_size * config.mesh_size * config.mesh_size;

    // Pre-calculate constants to avoid doing it inside the loop
    const double gamma_minus_1 = config.gamma - 1.0;
    const double inv_gamma_minus_1 = 1.0 / gamma_minus_1;

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
    double* metal = metal_density.array().data();

#pragma omp parallel for schedule(static)
    for (int i = 0; i < total_cells; ++i) {
        double local_rho = rho[i];

        if (local_rho > density_floor) {
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

            // DUAL ENERGY SWITCH
            // If internal energy is > x% of total energy, trust total energy.
            // A double has 53 bits of mantissa precision.
            // Guarantee at least half the bits survive the subtraction.
            constexpr double BIT_THRESHOLD = 0x1p-26;
            if (e_tot > 0.0 && (e_tot - ke) / e_tot > BIT_THRESHOLD) {
                ie_tracked[i] =
                    (e_tot - ke);  // Sync tracked IE to prevent drift
                p = gamma_minus_1 * ie_tracked[i];
            } else {
                // Hypersonic regime: Trust tracked IE.
                p = gamma_minus_1 * ie_tracked[i];
                en[i] = p * inv_gamma_minus_1 + ke;  // Sync Total Energy
            }

            if (p < pressure_floor) {
                p = pressure_floor;
                ie_tracked[i] = pressure_floor * inv_gamma_minus_1;
                en[i] = ie_tracked[i] + ke;
            }

            pr[i] = p;
            metal[i] = std::max(0.0, std::min(metal[i], local_rho));
        } else {
            // Vacuum cell handling
            vx[i] = 0.0;
            vy[i] = 0.0;
            vz[i] = 0.0;
            pr[i] = pressure_floor;
            ie_tracked[i] = pressure_floor * inv_gamma_minus_1;
            en[i] = ie_tracked[i];
            metal[i] = 0.0;
        }
    }
}

void GasGrid::compute_and_apply_fluxes(double dt) {
    double factor = dt / config.cell_size;

    for (int axis = 0; axis < 3; ++axis) {
        update_primitive_variables();
        solver.compute_fluxes(*this, axis, config.gamma);

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
        Grid3D dv(config.mesh_size);
        dv.data = v_n->roll(-1, axis).array() - v_n->roll(1, axis).array();

        int total_cells =
            config.mesh_size * config.mesh_size * config.mesh_size;

        double* d_rho = density.array().data();
        double* d_mn = m_n->array().data();
        double* d_mt1 = m_t1->array().data();
        double* d_mt2 = m_t2->array().data();
        double* d_en = energy.array().data();
        double* d_ie = internal_energy.array().data();
        double* d_metal = metal_density.array().data();

        double* p_arr = pressure.array().data();
        double* dv_arr = dv.array().data();

        const double* f_rho = solver.get_flux_density().array().data();
        const double* f_mn = solver.get_flux_mom_n().array().data();
        const double* f_mt1 = solver.get_flux_mom_t1().array().data();
        const double* f_mt2 = solver.get_flux_mom_t2().array().data();
        const double* f_en = solver.get_flux_energy().array().data();
        const double* f_ie = solver.get_flux_ie().array().data();
        const double* f_metal = solver.get_flux_metal().array().data();

        const double* f_rho_sh = solver.get_flux_density_sh().array().data();
        const double* f_mn_sh = solver.get_flux_mom_n_sh().array().data();
        const double* f_mt1_sh = solver.get_flux_mom_t1_sh().array().data();
        const double* f_mt2_sh = solver.get_flux_mom_t2_sh().array().data();
        const double* f_en_sh = solver.get_flux_energy_sh().array().data();
        const double* f_ie_sh = solver.get_flux_ie_sh().array().data();
        const double* f_metal_sh = solver.get_flux_metal_sh().array().data();

#pragma omp parallel for simd schedule(static)
        for (int i = 0; i < total_cells; ++i) {
            d_rho[i] -= factor * (f_rho[i] - f_rho_sh[i]);
            if (d_rho[i] < density_floor) d_rho[i] = density_floor;

            d_mn[i] -= factor * (f_mn[i] - f_mn_sh[i]);
            d_mt1[i] -= factor * (f_mt1[i] - f_mt1_sh[i]);
            d_mt2[i] -= factor * (f_mt2[i] - f_mt2_sh[i]);
            d_en[i] -= factor * (f_en[i] - f_en_sh[i]);

            d_ie[i] -= factor * (f_ie[i] - f_ie_sh[i]);
            d_ie[i] -= 0.5 * factor * p_arr[i] * dv_arr[i];  // PdV work

            d_metal[i] -= factor * (f_metal[i] - f_metal_sh[i]);
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
    Grid3D old_metal = metal_density;

    // RK2 Stage 1
    // Advance using the initial state: U^{(1)} = U^n + dt * L(U^n)
    compute_and_apply_fluxes(dt);

    // RK2 Stage 2
    // Advance again using the intermediate state: U^{(2)} = U^{(1)} + dt *
    // L(U^{(1)})
    compute_and_apply_fluxes(dt);

    // Final Averaging
    int total_cells = config.mesh_size * config.mesh_size * config.mesh_size;

    double* d_rho = density.array().data();
    double* d_mx = momentum_x.array().data();
    double* d_my = momentum_y.array().data();
    double* d_mz = momentum_z.array().data();
    double* d_en = energy.array().data();
    double* d_ie = internal_energy.array().data();
    double* d_metal = metal_density.array().data();

    double* o_rho = old_rho.array().data();
    double* o_mx = old_mx.array().data();
    double* o_my = old_my.array().data();
    double* o_mz = old_mz.array().data();
    double* o_en = old_en.array().data();
    double* o_ie = old_ie.array().data();
    double* o_metal = old_metal.array().data();

#pragma omp parallel for simd schedule(static)
    for (int i = 0; i < total_cells; ++i) {
        d_rho[i] = 0.5 * (o_rho[i] + d_rho[i]);
        d_mx[i] = 0.5 * (o_mx[i] + d_mx[i]);
        d_my[i] = 0.5 * (o_my[i] + d_my[i]);
        d_mz[i] = 0.5 * (o_mz[i] + d_mz[i]);
        d_en[i] = 0.5 * (o_en[i] + d_en[i]);
        d_ie[i] = 0.5 * (o_ie[i] + d_ie[i]);
        d_metal[i] = 0.5 * (o_metal[i] + d_metal[i]);
    }

    // Update primitive variables at the end so the gravity step has the
    // right values
    update_primitive_variables();
}

double GasGrid::get_cfl_timestep() const {
    if (config.hydro_method != HydroMethod::Eulerian)
        return std::numeric_limits<double>::infinity();

    Grid3D cs_sq(config.mesh_size);
    cs_sq.data = (config.gamma * pressure.array() / density.array());
    cs_sq.data = (density.array() > density_floor).select(cs_sq.array(), 0.0);
    cs_sq.data = (cs_sq.array() == cs_sq.array()).select(cs_sq.array(), 0.0);

    Grid3D v_mag(config.mesh_size);
    v_mag.data = (velocity_x.array().square() + velocity_y.array().square() +
                  velocity_z.array().square())
                     .sqrt();

    Grid3D signal_vel(config.mesh_size);
    signal_vel.data = v_mag.array() + cs_sq.array().sqrt();

    double max_signal_vel = signal_vel.maxCoeff();
    if (max_signal_vel < 1e-9) return std::numeric_limits<double>::infinity();

    return (config.cell_size / max_signal_vel) * config.hydro_courant_factor;
}

void GasGrid::apply_cooling(double dt, double a, Cooling& cooling) {
    if (!config.enable_cooling) return;

    double u_rad_floor = cooling.get_u_rad_floor(a, config);
    int total_cells = config.mesh_size * config.mesh_size * config.mesh_size;

    const double* d_rho = density.array().data();
    double* d_en = energy.array().data();
    double* d_ie = internal_energy.array().data();
    double* d_metal = metal_density.array().data();

    double total_radiated = 0.0;
    double total_photoheated = 0.0;
    size_t non_converged_count = 0;
    size_t total_cycles = 0;

#pragma omp parallel for schedule(static)                                 \
    reduction(+ : total_radiated, total_photoheated, non_converged_count, \
                  total_cycles)
    for (int i = 0; i < total_cells; ++i) {
        double local_rho = d_rho[i];

        if (local_rho > density_floor) {
            double local_Z_frac = d_metal[i] / local_rho;
            double u_current = d_ie[i] / local_rho;
            double u_initial = u_current;

            double t_evolved = 0.0;
            int cell_non_converged = 0;

            // Local Cell Subcycling
            while (t_evolved < dt) {
                // Determine a sub-step for this cell (e.g., max 10% energy
                // change)
                double du_dt = cooling.compute_du_dt(u_current, local_rho,
                                                     local_Z_frac, a, config);

                double dt_cell = (std::abs(du_dt) > 0.0)
                                     ? 0.1 * (u_current / std::abs(du_dt))
                                     : dt;

                // Don't overshoot the global hydro timestep
                dt_cell = std::min(dt_cell, dt - t_evolved);

                // Run the implicit solver for this tiny step
                int iters = 0;
                u_current = cooling.solve_cooling_implicit(
                    u_current, local_rho, local_Z_frac, a, dt_cell, u_rad_floor,
                    config, iters);

                if (iters >= Cooling::MAX_ITER) {
                    cell_non_converged++;
                }

                t_evolved += dt_cell;
                total_cycles++;
            }

            // Log if the cell struggled during any subcycle
            if (cell_non_converged > 0) {
                non_converged_count++;
            }

            // Calculate the total energy change over the step
            double delta_u = u_current - u_initial;

            if (std::abs(delta_u) > 0.0) {
                double delta_E_vol = delta_u * local_rho;
                // Add the delta directly
                d_ie[i] += delta_E_vol;
                d_en[i] += delta_E_vol;

                // Accumulate energy for logging
                if (delta_u < 0.0) {
                    // Net Cooling (Radiated away)
                    total_radiated -= delta_E_vol * config.cell_volume;
                } else {
                    // Net Heating (Injected by UVB)
                    total_photoheated += delta_E_vol * config.cell_volume;
                }
            }
        }
    }

    this->cooling_total_cycles = total_cycles / total_cells;
    this->cooling_failed_cells = non_converged_count;
    this->accumulated_radiated_energy += total_radiated;
    this->accumulated_photoheating_energy += total_photoheated;

    update_primitive_variables();
}

double GasGrid::get_cooling_timestep(double a, Cooling& cooling) const {
    if (!config.enable_cooling) return std::numeric_limits<double>::infinity();

    double min_dt_cool = std::numeric_limits<double>::infinity();
    int total_cells = config.mesh_size * config.mesh_size * config.mesh_size;

    const double* d_rho = density.array().data();
    const double* d_ie = internal_energy.array().data();
    const double* d_metal = metal_density.array().data();

    // Calculate the physical cooling floor
    double u_rad_floor = cooling.get_u_rad_floor(a, config);

#pragma omp parallel for reduction(min : min_dt_cool)
    for (int i = 0; i < total_cells; ++i) {
        if (d_rho[i] > density_floor) {
            double u = d_ie[i] / d_rho[i];

            if (u <= u_rad_floor) continue;

            double Z_frac = d_metal[i] / d_rho[i];
            double du_dt =
                cooling.compute_du_dt(u, d_rho[i], Z_frac, a, config);
            if (std::abs(du_dt) > 0.0) {
                // Restrict timestep so internal energy changes at most 10%
                double dt_cool = 0.1 * (u / std::abs(du_dt));
                if (dt_cool < min_dt_cool) {
                    min_dt_cool = dt_cool;
                }
            }
        }
    }

    return min_dt_cool;
}

Grid3D GasGrid::compute_thermal_timescale(double a,
                                          const Cooling& cooling) const {
    Grid3D t_therm_grid(config.mesh_size);

    double* t_therm = t_therm_grid.array().data();
    const double* d_rho = density.array().data();
    const double* d_ie = internal_energy.array().data();
    const double* d_metal = metal_density.array().data();

    const int total_cells =
        config.mesh_size * config.mesh_size * config.mesh_size;

#pragma omp parallel for schedule(static)
    for (int i = 0; i < total_cells; ++i) {
        if (d_rho[i] > density_floor) {
            double u = d_ie[i] / d_rho[i];
            double Z_frac = d_metal[i] / d_rho[i];

            // Calculate the instantaneous change in internal energy
            double du_dt =
                cooling.compute_du_dt(u, d_rho[i], Z_frac, a, config);

            // Prevent division by zero if the cell is in equilibrium
            if (std::abs(du_dt) > 1e-30) {
                t_therm[i] = u / du_dt;
            } else {
                t_therm[i] = std::numeric_limits<double>::infinity();
            }
        } else {
            // Vacuum cells have effectively infinite thermal timescales
            t_therm[i] = std::numeric_limits<double>::infinity();
        }
    }

    return t_therm_grid;
}
