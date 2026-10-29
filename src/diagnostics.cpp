#include "diagnostics.h"

#include "mfm.h"
#include "state.h"

ScopedTimer::ScopedTimer(Diagnostics& d, TimerRegion r) : diag(d), region(r) {
    start_time = std::chrono::high_resolution_clock::now();
}
ScopedTimer::~ScopedTimer() {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    diag.add_time(region, elapsed.count());
}

void Diagnostics::add_prof_time(ProfRegion region, double time_sec) {
    accumulated_prof[static_cast<size_t>(region)] += time_sec;
}

void Diagnostics::add_time(TimerRegion region, double time_sec) {
    accumulated_times[static_cast<size_t>(region)] += time_sec;
}

void Diagnostics::increment_cycle() {
    cycle++;
    accumulated_cycles++;
}

double Diagnostics::get_average(TimerRegion region) const {
    if (accumulated_cycles == 0) return 0.0;
    return accumulated_times[static_cast<size_t>(region)] / accumulated_cycles;
}

double Diagnostics::get_io_time() const {
    // I/O is usually a single spike, so we just return the accumulated time
    return accumulated_times[static_cast<size_t>(TimerRegion::IO)];
}

double Diagnostics::get_average_overhead() const {
    double total_physics =
        get_average(TimerRegion::PM) + get_average(TimerRegion::PP) +
        get_average(TimerRegion::Hydro) + get_average(TimerRegion::Cool) +
        get_average(TimerRegion::Tree) + get_average(TimerRegion::Density);
    return get_average(TimerRegion::Step) - total_physics;
}

double Diagnostics::get_prof_average(ProfRegion region) const {
    if (accumulated_cycles == 0) return 0.0;
    return accumulated_prof[static_cast<size_t>(region)] / accumulated_cycles;
}

void Diagnostics::reset_accumulators() {
    accumulated_times.fill(0.0);
    accumulated_prof.fill(0.0);
    accumulated_cycles = 0;
}

struct HydroConservationError {
    Eigen::Vector3d momentum_error_rate;  // Net un-conserved force
    double energy_error_rate;             // Net un-conserved power
};

static HydroConservationError diagnose_hydro_conservation_error(
    const GasParticleSystem& gas) {
    double net_fx = 0.0;
    double net_fy = 0.0;
    double net_fz = 0.0;
    double net_edot = 0.0;

    // Summing m * a_hydro gives the aggregate of all Riemann face fluxes
#pragma omp parallel for reduction(+ : net_fx, net_fy, net_fz, net_edot) \
    schedule(static)
    for (size_t i = 0; i < gas.num_particles; ++i) {
        net_fx += gas.mass[i] * gas.hydro_acc_x[i];
        net_fy += gas.mass[i] * gas.hydro_acc_y[i];
        net_fz += gas.mass[i] * gas.hydro_acc_z[i];
        net_edot += gas.mass[i] * gas.de_dt[i];
    }

    HydroConservationError error;
    error.momentum_error_rate = Eigen::Vector3d(net_fx, net_fy, net_fz);
    error.energy_error_rate = net_edot;

    return error;
}

void Diagnostics::compute_dm_metrics(const ParticleSystem& dm) {
    double sum_mass = 0.0, sum_px = 0.0, sum_py = 0.0, sum_pz = 0.0;
    double k_dm_code = 0.0;

#pragma omp parallel for reduction(+ : sum_mass, sum_px, sum_py, sum_pz, \
                                       k_dm_code)
    for (size_t i = 0; i < dm.num_particles; ++i) {
        double m = dm.mass[i];
        double vx = dm.vel_x[i];
        double vy = dm.vel_y[i];
        double vz = dm.vel_z[i];

        sum_mass += m;
        sum_px += m * vx;
        sum_py += m * vy;
        sum_pz += m * vz;
        k_dm_code += 0.5 * m * (vx * vx + vy * vy + vz * vz);
    }

    this->total_mass_dm = sum_mass;
    this->total_momentum_dm = {sum_px, sum_py, sum_pz};
    this->ke_dm = k_dm_code;

    // Evaluate DM Energy Error
    if (!dm_energy_initialized) {
        initial_dm_energy = ke_dm;
        dm_energy_initialized = true;
    }

    double delta_ke_dm = ke_dm - initial_dm_energy;
    double dm_absolute_error = delta_ke_dm - dm.accumulated_gravitational_work +
                               dm.accumulated_expansion_work;
    double denominator =
        (initial_dm_energy > 1e-12) ? initial_dm_energy : ke_dm;

    this->dm_energy_err = (std::abs(denominator) > 1e-12)
                              ? std::abs(dm_absolute_error / denominator)
                              : 0.0;
}

double Diagnostics::calculate_hydro_energy_error(double current_energy,
                                                 double w_grav, double w_exp,
                                                 double e_rad, double e_heat) {
    if (!energy_initialized) {
        initial_gas_energy = current_energy;
        energy_initialized = true;
    }

    double delta_e = current_energy - initial_gas_energy;
    double absolute_error = delta_e - w_grav + w_exp + e_rad - e_heat;
    double denominator = std::abs(initial_gas_energy);

    return (denominator > 1e-12) ? absolute_error / denominator
                                 : absolute_error;
}

void Diagnostics::compute_eulerian_metrics(const GasGrid& gas,
                                           const Config& config) {
    this->max_gas_density = gas.get_density().maxCoeff();
    this->max_gas_pressure = gas.get_pressure().maxCoeff();
    this->max_gas_velocity = (gas.get_velocity_x().array().square() +
                              gas.get_velocity_y().array().square() +
                              gas.get_velocity_z().array().square())
                                 .sqrt()
                                 .maxCoeff();

    this->non_converged_cooling_cells = gas.get_cooling_failed_cells();
    this->total_mass_gas = gas.get_density().sum() * config.cell_volume;
    this->total_momentum_gas = {
        gas.get_momentum_x().sum() * config.cell_volume,
        gas.get_momentum_y().sum() * config.cell_volume,
        gas.get_momentum_z().sum() * config.cell_volume};

    int total_cells = gas.get_density().size();
    Grid3D ke_gas_density(config.mesh_size);
    ke_gas_density.data = 0.5 * (gas.get_momentum_x().array().square() +
                                 gas.get_momentum_y().array().square() +
                                 gas.get_momentum_z().array().square());

#pragma omp parallel for simd schedule(static)
    for (int i = 0; i < total_cells; ++i) {
        ke_gas_density.data[i] /= std::max(1e-12, gas.get_density().data[i]);
    }

    this->ke_gas = ke_gas_density.sum() * config.cell_volume;
    this->ie_gas = (gas.get_pressure().array() / (config.gamma - 1.0)).sum() *
                   config.cell_volume;

    this->total_radiated_energy = gas.get_accumulated_radiated_energy();
    this->total_heated_energy = gas.get_accumulated_photoheating_energy();

    this->energy_err = calculate_hydro_energy_error(
        this->ke_gas + this->ie_gas, gas.get_accumulated_gravitational_work(),
        gas.get_accumulated_expansion_work(), this->total_radiated_energy,
        this->total_heated_energy);
}

void Diagnostics::compute_mfm_metrics(GasParticleSystem& gas,
                                      const Config& config) {
    double max_rho = 0.0, max_p = 0.0, max_v = 0.0;
    double sum_mass = 0.0, sum_px = 0.0, sum_py = 0.0, sum_pz = 0.0;
    double sum_ke = 0.0, sum_ie = 0.0, total_sim_volume = 0.0;

#pragma omp parallel for reduction(max : max_rho, max_p, max_v)     \
    reduction(+ : sum_mass, sum_px, sum_py, sum_pz, sum_ke, sum_ie, \
                  total_sim_volume)
    for (size_t i = 0; i < gas.num_particles; ++i) {
        max_rho = std::max(max_rho, gas.rho[i]);
        max_p = std::max(max_p, gas.pressure[i]);

        double v_sq = gas.vel_x[i] * gas.vel_x[i] +
                      gas.vel_y[i] * gas.vel_y[i] + gas.vel_z[i] * gas.vel_z[i];
        max_v = std::max(max_v, std::sqrt(v_sq));

        double m = gas.mass[i];
        sum_mass += m;
        sum_px += m * gas.vel_x[i];
        sum_py += m * gas.vel_y[i];
        sum_pz += m * gas.vel_z[i];

        sum_ke += 0.5 * m * v_sq;
        sum_ie += m * gas.u[i];
        total_sim_volume += m / std::max(gas.rho[i], 1e-12);
    }

    this->max_gas_density = max_rho;
    this->max_gas_pressure = max_p;
    this->max_gas_velocity = max_v;
    this->total_mass_gas = sum_mass;
    this->total_momentum_gas = {sum_px, sum_py, sum_pz};
    this->ke_gas = sum_ke;
    this->ie_gas = sum_ie;
    this->non_converged_cooling_cells = gas.cooling_failed_cells;

    // Extracted Energy Error
    this->total_radiated_energy = gas.accumulated_radiated_energy;
    this->total_heated_energy = gas.accumulated_photoheating_energy;
    this->energy_err = calculate_hydro_energy_error(
        sum_ke + sum_ie, gas.accumulated_gravitational_work,
        gas.accumulated_expansion_work, gas.accumulated_radiated_energy,
        gas.accumulated_photoheating_energy);

    // MFM Specific Diagnostics
    this->MFM_total_particles = gas.num_particles;
    this->MFM_h_non_converged_cases = gas.non_converged_h_cases;
    this->MFM_not_enough_nbg_cases = gas.not_enough_neighbors_cases;
    this->MFM_sph_fallback_cases = gas.sph_fallback_cases;
    this->MFM_riemann_failed_cases = gas.riemann_failed_cases;

    this->percent_particles_updated = (config.individual_particle_timesteps &&
                                       gas.active_particles_num_cycles > 0)
                                          ? 100.0 *
                                                gas.active_particles_fraction /
                                                gas.active_particles_num_cycles
                                          : 100.0;

    double cycle_div = std::max(1, cycle);
    this->MFM_ill_percent =
        100.0 * (static_cast<double>(gas.ill_conditioned_cases) / cycle_div) /
        MFM_total_particles;
    this->MFM_clamp_percent =
        100.0 * (static_cast<double>(gas.clamped_h_cases) / cycle_div) /
        MFM_total_particles;
    this->MFM_inc_percent = 100.0 *
                            (static_cast<double>(gas.neighbor_increased_cases) /
                             std::max(1, accumulated_cycles)) /
                            MFM_total_particles;

    double box_volume =
        config.domain_size.x * config.domain_size.y * config.domain_size.z;
    this->volume_error = (total_sim_volume / box_volume) - 1.0;

    auto hce = diagnose_hydro_conservation_error(gas);
    this->MFM_flux_energy_error = hce.energy_error_rate;
    this->MFM_flux_momentum_error = hce.momentum_error_rate.norm();

    // Reset counters
    gas.active_particles_fraction = 0.0;
    gas.active_particles_num_cycles = 0;
    gas.ill_conditioned_cases = 0;
    gas.clamped_h_cases = 0;
    gas.neighbor_increased_cases = 0;
}

double Diagnostics::compute_potential_energy(const Grid3D& total_rho,
                                             const Grid3D& phi,
                                             const Config& config) {
    double w_code = 0.0;
    int total_cells = total_rho.data.size();
    double n3_inv = 1.0 / static_cast<double>(total_cells);

#pragma omp parallel for reduction(+ : w_code)
    for (int i = 0; i < total_cells; ++i) {
        w_code += 0.5 * total_rho.data[i] * (phi.data[i] * n3_inv) *
                  config.cell_volume;
    }
    return w_code;
}

void Diagnostics::update_physics(const SimState& state, const TimestepInfo& ts,
                                 const Config& config) {
    // Time and Timestep Metrics
    this->sim_time = state.total_time;
    this->scale_factor = state.scale_factor;
    this->dt_cfl = ts.dt_hydro;
    this->dt_gravity = ts.dt_grav;
    this->dt_cool = ts.dt_cool;
    this->dt_final = ts.dt_macro;

    // Compute Independent Physics Domains
    compute_dm_metrics(state.dm);
    this->pe_total =
        compute_potential_energy(state.pm.rho, state.pm.phi, config);

    // Compute Solver-Specific Hydro Metrics
    if (config.hydro_method == HydroMethod::Eulerian) {
        compute_eulerian_metrics(*state.gas, config);
    } else if (config.hydro_method == HydroMethod::MFM) {
        compute_mfm_metrics(*state.mfm_gas, config);
    }

    // Global System Error Metrics
    this->total_momentum = (total_momentum_dm + total_momentum_gas).norm();

    if (initial_mass == 0.0) {
        initial_mass = total_mass_dm + total_mass_gas;
    }
    this->mass_err = total_mass_dm + total_mass_gas - initial_mass;
}
