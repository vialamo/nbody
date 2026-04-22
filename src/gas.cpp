#include "gas.h"
#include <omp.h>

RiemannSolver::RiemannSolver(int mesh_size)
    : density(mesh_size),
      mom_n(mesh_size),
      mom_t1(mesh_size),
      mom_t2(mesh_size),
      energy(mesh_size),
      v_n(mesh_size),
      v_t1(mesh_size),
      v_t2(mesh_size),
      pressure(mesh_size),
      rho_L(mesh_size),
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
      S_R_minus_S_L(mesh_size),
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
      flux_energy(mesh_size) {}

Grid3D RiemannSolver::solve_hll(const Grid3D& FL, const Grid3D& FR,
                                const Grid3D& UL, const Grid3D& UR) {
    Grid3D flux(density.n);
    flux.data =
        (S_L.array() >= 0)
            .select(FL.array(),
                    (S_R.array() <= 0)
                        .select(FR.array(), (S_R.array() * FL.array() -
                                             S_L.array() * FR.array() +
                                             S_L.array() * S_R.array() *
                                                 (UR.array() - UL.array())) /
                                                S_R_minus_S_L.array()));
    return flux;
}

void RiemannSolver::compute_fluxes(const GasGrid& grid, int axis,
                                   double gamma) {
    // --- 1. Map Global Grid to Local Coordinate System (Normal/Transverse) ---
    this->density = grid.get_density();
    this->energy = grid.get_energy();
    this->pressure = grid.get_pressure();

    if (axis == 0) {  // X-sweep: Normal is X, Transverse are Y and Z
        mom_n = grid.get_momentum_x();
        v_n = grid.get_velocity_x();
        mom_t1 = grid.get_momentum_y();
        v_t1 = grid.get_velocity_y();
        mom_t2 = grid.get_momentum_z();
        v_t2 = grid.get_velocity_z();
    } else if (axis == 1) {  // Y-sweep: Normal is Y, Transverse are Z and X
        mom_n = grid.get_momentum_y();
        v_n = grid.get_velocity_y();
        mom_t1 = grid.get_momentum_z();
        v_t1 = grid.get_velocity_z();
        mom_t2 = grid.get_momentum_x();
        v_t2 = grid.get_velocity_x();
    } else {  // Z-sweep: Normal is Z, Transverse are X and Y
        mom_n = grid.get_momentum_z();
        v_n = grid.get_velocity_z();
        mom_t1 = grid.get_momentum_x();
        v_t1 = grid.get_velocity_x();
        mom_t2 = grid.get_momentum_y();
        v_t2 = grid.get_velocity_y();
    }

    // --- 2. Reconstruct Left (L) and Right (R) States ---
    // In this first-order scheme, L is the current cell and R is the neighbor
    // in +axis direction
    const int rollDir = -1;
    rho_L = density;
    rho_R = rho_L.roll(rollDir, axis);
    p_L = pressure;
    p_R = p_L.roll(rollDir, axis);
    vn_L = v_n;
    vn_R = vn_L.roll(rollDir, axis);
    vt1_L = v_t1;
    vt1_R = vt1_L.roll(rollDir, axis);
    vt2_L = v_t2;
    vt2_R = vt2_L.roll(rollDir, axis);
    E_L = energy;
    E_R = E_L.roll(rollDir, axis);

    mom_n_L = mom_n;
    mom_n_R = mom_n_L.roll(rollDir, axis);
    mom_t1_L = mom_t1;
    mom_t1_R = mom_t1_L.roll(rollDir, axis);
    mom_t2_L = mom_t2;
    mom_t2_R = mom_t2_L.roll(rollDir, axis);

    // --- 3. Compute Sound Speeds and Wave Signal Speeds ---
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

    // --- 4. Evaluate Physical Fluxes for Left and Right states ---
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

    // Preparation for division in HLL formula
    S_R_minus_S_L.data = S_R.array() - S_L.array();
    S_R_minus_S_L.data =
        (S_R_minus_S_L.array().abs() < 1e-9).select(1e-9, S_R_minus_S_L.data);

    // --- 5. Solve for Intercell Fluxes using HLL ---
    flux_density = solve_hll(F_dens_L, F_dens_R, rho_L, rho_R);
    flux_mom_n = solve_hll(F_momn_L, F_momn_R, mom_n_L, mom_n_R);
    flux_mom_t1 = solve_hll(F_momt1_L, F_momt1_R, mom_t1_L, mom_t1_R);
    flux_mom_t2 = solve_hll(F_momt2_L, F_momt2_R, mom_t2_L, mom_t2_R);
    flux_energy = solve_hll(F_en_L, F_en_R, E_L, E_R);

    // --- 6. Pre-calculate shifted fluxes for the grid update ---
    flux_density_sh = flux_density.roll(1, axis);
    flux_mom_n_sh = flux_mom_n.roll(1, axis);
    flux_mom_t1_sh = flux_mom_t1.roll(1, axis);
    flux_mom_t2_sh = flux_mom_t2.roll(1, axis);
    flux_energy_sh = flux_energy.roll(1, axis);
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
      solver(conf.MESH_SIZE),
      config(conf) {}

void GasGrid::update_primitive_variables() {
    int total_cells = config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;

    #pragma omp parallel for
    for (int i = 0; i < total_cells; ++i) {
        if (density.data[i] > 1e-12) {
            velocity_x.data[i] = momentum_x.data[i] / density.data[i];
            velocity_y.data[i] = momentum_y.data[i] / density.data[i];
            velocity_z.data[i] = momentum_z.data[i] / density.data[i];
        } else {
            velocity_x.data[i] = velocity_y.data[i] = velocity_z.data[i] = 0.0;
        }
    }

    Grid3D kin_energy(config.MESH_SIZE);
    kin_energy.data =
        0.5 * (momentum_x.array().square() + momentum_y.array().square() +
               momentum_z.array().square());

    #pragma omp parallel for
    for (int i = 0; i < total_cells; ++i) {
        if (density.data[i] > 1e-12)
            kin_energy.data[i] /= density.data[i];
        else
            kin_energy.data[i] = 0.0;
    }

    pressure.data = (config.GAMMA - 1.0) * (energy.data - kin_energy.data);
    pressure.data = (pressure.array() < 1e-12).select(1e-12, pressure.data);
    energy.data = pressure.data / (config.GAMMA - 1.0) + kin_energy.data;
}

void GasGrid::hydro_step(double dt) {
    double factor = dt / config.CELL_SIZE;

    for (int axis = 0; axis < 3; ++axis) {
        update_primitive_variables();
        solver.compute_fluxes(*this, axis, config.GAMMA);

        // Helper references to correctly update original coordinates
        Grid3D* m_n;
        Grid3D* m_t1;
        Grid3D* m_t2;
        if (axis == 0) {
            m_n = &momentum_x;
            m_t1 = &momentum_y;
            m_t2 = &momentum_z;
        } else if (axis == 1) {
            m_n = &momentum_y;
            m_t1 = &momentum_z;
            m_t2 = &momentum_x;
        } else {
            m_n = &momentum_z;
            m_t1 = &momentum_x;
            m_t2 = &momentum_y;
        }

        density.array() -= factor * (solver.get_flux_density().array() -
                                     solver.get_flux_density_sh().array());

        // Enforce the density floor to prevent NaN crashes in voids
        density.data = (density.array() < 1e-12).select(1e-12, density.data);

        m_n->array() -= factor * (solver.get_flux_mom_n().array() -
                                  solver.get_flux_mom_n_sh().array());

        m_t1->array() -= factor * (solver.get_flux_mom_t1().array() -
                                   solver.get_flux_mom_t1_sh().array());

        m_t2->array() -= factor * (solver.get_flux_mom_t2().array() -
                                   solver.get_flux_mom_t2_sh().array());

        energy.array() -= factor * (solver.get_flux_energy().array() -
                                    solver.get_flux_energy_sh().array());
    }
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