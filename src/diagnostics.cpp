#include "diagnostics.h"

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
        get_average(TimerRegion::Hydro) + get_average(TimerRegion::Cool);
    return get_average(TimerRegion::Step) - total_physics;
}

double Diagnostics::get_prof_average(ProfRegion region) const {
    if (accumulated_cycles == 0) return 0.0;
    return accumulated_prof[static_cast<size_t>(region)] / accumulated_cycles;
}

void Diagnostics::reset_accumulators() {
    accumulated_times.fill(0.0);
    accumulated_substeps.fill(0);
    accumulated_prof.fill(0.0);
    accumulated_cycles = 0;
}

void Diagnostics::update_physics(const SimState& state, const TimestepInfo& ts,
                                 const Config& config) {
    this->sim_time = state.total_time;
    this->scale_factor = state.scale_factor;
    this->dt_cfl = ts.dt_hydro;
    this->dt_gravity = ts.dt_grav;
    this->dt_cool = ts.dt_cool;
    this->dt_final = ts.dt_macro;

    double k_gas_code = 0.0, u_gas_code = 0.0;

    if (config.use_hydro) {
        this->max_gas_density = state.gas.get_density().maxCoeff();
        this->max_gas_pressure = state.gas.get_pressure().maxCoeff();
        this->max_gas_velocity = (state.gas.get_velocity_x().array().square() +
                                  state.gas.get_velocity_y().array().square() +
                                  state.gas.get_velocity_z().array().square())
                                     .sqrt()
                                     .maxCoeff();
        this->non_converged_cooling_cells =
            state.gas.get_cooling_failed_cells();

        this->total_mass_gas =
            state.gas.get_density().sum() * config.cell_volume;
        this->total_momentum_gas.x() =
            state.gas.get_momentum_x().sum() * config.cell_volume;
        this->total_momentum_gas.y() =
            state.gas.get_momentum_y().sum() * config.cell_volume;
        this->total_momentum_gas.z() =
            state.gas.get_momentum_z().sum() * config.cell_volume;

        int total_cells =
            config.mesh_size * config.mesh_size * config.mesh_size;
        Grid3D ke_gas_density(config.mesh_size);
        ke_gas_density.data =
            0.5 * (state.gas.get_momentum_x().array().square() +
                   state.gas.get_momentum_y().array().square() +
                   state.gas.get_momentum_z().array().square());

#pragma omp parallel for simd schedule(static)
        for (int i = 0; i < total_cells; ++i) {
            ke_gas_density.data[i] /=
                std::max(1e-12, state.gas.get_density().data[i]);
        }
        k_gas_code = ke_gas_density.sum() * config.cell_volume;
        u_gas_code =
            (state.gas.get_pressure().array() / (config.gamma - 1.0)).sum() *
            config.cell_volume;
    }

    const auto& dm = state.dm;
    double sum_mass = 0.0, sum_px = 0.0, sum_py = 0.0, sum_pz = 0.0,
           k_dm_code = 0.0;

#pragma omp parallel for reduction(+ : sum_mass, sum_px, sum_py, sum_pz, \
                                       k_dm_code)
    for (size_t i = 0; i < dm.num_particles; ++i) {
        sum_mass += dm.mass[i];
        sum_px += dm.mass[i] * dm.vel_x[i];
        sum_py += dm.mass[i] * dm.vel_y[i];
        sum_pz += dm.mass[i] * dm.vel_z[i];

        k_dm_code += 0.5 * dm.mass[i] *
                     (dm.vel_x[i] * dm.vel_x[i] + dm.vel_y[i] * dm.vel_y[i] +
                      dm.vel_z[i] * dm.vel_z[i]);
    }

    this->total_mass_dm = sum_mass;
    this->total_momentum_dm = {sum_px, sum_py, sum_pz};
    this->total_momentum = (total_momentum_dm + total_momentum_gas).norm();
    this->mass_err = total_mass_dm + total_mass_gas - config.total_mass;

    // Grid PM potential energy
    double w_code = 0.0;
    if (config.enable_energy_diagnostics) {
        int total_cells =
            config.mesh_size * config.mesh_size * config.mesh_size;
        double n3_inv = 1.0 / static_cast<double>(total_cells);
#pragma omp parallel for reduction(+ : w_code)
        for (int i = 0; i < total_cells; ++i) {
            w_code += 0.5 * state.total_rho.data[i] *
                      (state.phi.data[i] * n3_inv) * config.cell_volume;
        }
    }

    // Raw code units
    this->ke_dm = k_dm_code;
    this->ke_gas = k_gas_code;
    this->ie_gas = u_gas_code;
    this->pe_total = w_code;

    // Evaluate Gas Error
    if (config.use_hydro) {
        double current_gas_energy = k_gas_code + u_gas_code;

        if (!energy_initialized) {
            initial_gas_energy = current_gas_energy;
            energy_initialized = true;
        }

        double delta_e = current_gas_energy - initial_gas_energy;
        double w_grav = state.gas.get_accumulated_gravitational_work();
        double w_exp = state.gas.get_accumulated_expansion_work();
        double e_rad = state.gas.get_accumulated_radiated_energy();
        double e_heat = state.gas.get_accumulated_photoheating_energy();

        double absolute_error = delta_e - w_grav + w_exp + e_rad - e_heat;
        double denominator = std::abs(initial_gas_energy);
        this->energy_err =
            (denominator > 1e-12) ? absolute_error / denominator : 0.0;

        total_radiated_energy = e_rad;
        total_heated_energy = e_heat;
    } else {
        this->energy_err = 0.0;
    }
}