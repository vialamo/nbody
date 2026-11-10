#include "diagnostics.h"

#include "gas.h"
#include "state.h"

namespace {
struct HydroConservationError {
    Eigen::Vector3d momentum_error_rate;  // Net un-conserved force
    double energy_error_rate;             // Net un-conserved power
};

HydroConservationError diagnose_hydro_conservation_error(
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
}  // namespace

ScopedTimer::ScopedTimer(Diagnostics& d, TimerRegion r) : diag(d), region(r) {
    start_time = std::chrono::high_resolution_clock::now();
}
ScopedTimer::~ScopedTimer() {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    diag.add_time(region, elapsed.count());
}

void Diagnostics::add_prof_time(GPUProfRegion region, double time_sec) {
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

double Diagnostics::get_prof_average(GPUProfRegion region) const {
    if (accumulated_cycles == 0) return 0.0;
    return accumulated_prof[static_cast<size_t>(region)] / accumulated_cycles;
}

void Diagnostics::reset_accumulators() {
    accumulated_times.fill(0.0);
    accumulated_prof.fill(0.0);
    accumulated_cycles = 0;
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
    double dm_absolute_error =
        delta_ke_dm - dm.budget.gravitational_work + dm.budget.expansion_work;
    double denominator =
        (initial_dm_energy > 1e-12) ? initial_dm_energy : ke_dm;

    this->dm_energy_err = (std::abs(denominator) > 1e-12)
                              ? std::abs(dm_absolute_error / denominator)
                              : 0.0;
}

double Diagnostics::calculate_hydro_energy_error(double current_energy,
                                                 const EnergyBudget& budget) {
    if (!gas_energy_initialized) {
        initial_gas_energy = current_energy;
        gas_energy_initialized = true;
    }

    double delta_e = current_energy - initial_gas_energy;
    double absolute_error = delta_e - budget.gravitational_work +
                            budget.expansion_work + budget.radiated -
                            budget.photoheated;
    double denominator = std::abs(initial_gas_energy);

    return (denominator > 1e-12) ? absolute_error / denominator
                                 : absolute_error;
}

void Diagnostics::compute_gas_metrics(const GasParticleSystem& gas,
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
    this->non_converged_cooling_particles = gas.diag.cooling_failed;

    // Extracted Energy Error
    this->total_radiated_energy = gas.budget.radiated;
    this->total_heated_energy = gas.budget.photoheated;
    this->gas_energy_err =
        calculate_hydro_energy_error(sum_ke + sum_ie, gas.budget);

    // MFM Specific Diagnostics
    const size_t num_gas_particles = gas.num_particles;
    this->MFM_h_non_converged_cases = gas.diag.non_converged_h;
    this->MFM_not_enough_nbg_cases = gas.diag.not_enough_neighbors;
    this->MFM_sph_fallback_cases = gas.diag.sph_fallback;
    this->MFM_riemann_failed_cases = gas.diag.riemann_failed;

    this->percent_particles_updated =
        (config.individual_particle_timesteps &&
         gas.diag.active_particles_num_cycles > 0)
            ? 100.0 * gas.diag.active_particles_fraction /
                  gas.diag.active_particles_num_cycles
            : 100.0;

    double cycle_div = std::max(1, cycle);
    if (num_gas_particles > 0) {
        this->MFM_ill_percent =
            100.0 *
            (static_cast<double>(gas.diag.ill_conditioned) / cycle_div) /
            num_gas_particles;
        this->MFM_clamp_percent =
            100.0 * (static_cast<double>(gas.diag.clamped_h) / cycle_div) /
            num_gas_particles;
        this->MFM_inc_percent =
            100.0 *
            (static_cast<double>(gas.diag.neighbor_increased) /
             std::max(1, accumulated_cycles)) /
            num_gas_particles;
    }

    double box_volume =
        config.domain_size.x * config.domain_size.y * config.domain_size.z;
    this->volume_error = (total_sim_volume / box_volume) - 1.0;

    auto hce = diagnose_hydro_conservation_error(gas);
    this->MFM_flux_energy_error = hce.energy_error_rate;
    this->MFM_flux_momentum_error = hce.momentum_error_rate.norm();
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
    this->softening_length = state.gravity.softening_length;

    // Compute Independent Physics Domains
    compute_dm_metrics(state.dm);
    this->pe_total =
        compute_potential_energy(state.pm.rho, state.pm.phi, config);

    // Compute Solver-Specific Hydro Metrics
    compute_gas_metrics(state.gas, config);

    // Global System Error Metrics
    this->total_momentum = (total_momentum_dm + total_momentum_gas).norm();

    if (initial_mass == 0.0) {
        initial_mass = total_mass_dm + total_mass_gas;
    }
    this->mass_err = total_mass_dm + total_mass_gas - initial_mass;
}
