#include "diagnostics.h"

ScopedTimer::ScopedTimer(Diagnostics& d, TimerRegion r) : diag(d), region(r) {
    start_time = std::chrono::high_resolution_clock::now();
}
ScopedTimer::~ScopedTimer() {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    diag.add_time(region, elapsed.count());
}

void Diagnostics::add_time(TimerRegion region, double time_sec) {
    accumulated_times[static_cast<size_t>(region)] += time_sec;
}

void Diagnostics::increment_cycle() { accumulated_cycles++; }

double Diagnostics::get_average(TimerRegion region) const {
    if (accumulated_cycles == 0) return 0.0;
    return accumulated_times[static_cast<size_t>(region)] / accumulated_cycles;
}

double Diagnostics::get_io_time() const {
    // I/O is usually a single spike, so we just return the raw accumulated time
    return accumulated_times[static_cast<size_t>(TimerRegion::IO)];
}

double Diagnostics::get_average_overhead() const {
    double total_physics = get_average(TimerRegion::PM) +
                           get_average(TimerRegion::PP) +
                           get_average(TimerRegion::Hydro);
    return get_average(TimerRegion::Step) - total_physics;
}

void Diagnostics::reset_accumulators() {
    accumulated_times.fill(0.0);
    accumulated_cycles = 0;
}

double Diagnostics::total_mass() const {
    return total_mass_dm + total_mass_gas;
}
double Diagnostics::total_energy() const {
    return ke_dm + ke_gas + ie_gas + pe_total;
}

void Diagnostics::update_physics(const SimState& state, double dt,
                                 const Config& config) {
    // Basic State
    this->sim_time = state.total_time;
    this->scale_factor = state.scale_factor;

    // Stability
    this->dt_cfl = state.gas.get_cfl_timestep();
    this->dt_gravity = state.dm.get_gravity_timestep(config);
    this->dt_final = dt;

    if (config.USE_HYDRO) {
        this->max_gas_density = state.gas.get_density().maxCoeff();
        this->max_gas_pressure = state.gas.get_pressure().maxCoeff();
        this->max_gas_velocity = (state.gas.get_velocity_x().array().square() +
                                  state.gas.get_velocity_y().array().square() +
                                  state.gas.get_velocity_z().array().square())
                                     .sqrt()
                                     .maxCoeff();
    }

    // Conservation (Particles)
    const auto& dm = state.dm;
    double sum_mass = 0.0, sum_px = 0.0, sum_py = 0.0, sum_pz = 0.0;

#pragma omp parallel for reduction(+ : sum_mass, sum_px, sum_py, sum_pz)
    for (size_t i = 0; i < dm.num_particles; ++i) {
        sum_mass += dm.mass[i];
        sum_px += dm.mass[i] * dm.vel_x[i];
        sum_py += dm.mass[i] * dm.vel_y[i];
        sum_pz += dm.mass[i] * dm.vel_z[i];
    }

    // Reset before adding
    this->total_mass_dm = sum_mass;
    this->total_momentum_dm = {sum_px, sum_py, sum_pz};
    this->ke_dm = 0.0;
    this->pe_total = 0.0;

    if (config.ENABLE_ENERGY_DIAGNOSTICS) {
        this->ke_dm = state.dm.calculate_kinetic_energy(state.scale_factor);

        double system_pe = 0.0;
        int total_cells =
            config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;

#pragma omp parallel for reduction(+ : system_pe)
        for (int i = 0; i < total_cells; ++i) {
            system_pe += 0.5 * state.total_rho.data[i] * state.phi.data[i] *
                         config.CELL_VOLUME;
        }
        this->pe_total = system_pe / state.scale_factor;
    }

    // Conservation (Gas)
    if (config.USE_HYDRO) {
        this->total_mass_gas =
            state.gas.get_density().sum() * config.CELL_VOLUME;
        this->total_momentum_gas.x() = state.gas.get_momentum_x().sum();
        this->total_momentum_gas.y() = state.gas.get_momentum_y().sum();
        this->total_momentum_gas.z() = state.gas.get_momentum_z().sum();

        Grid3D ke_gas_density(config.MESH_SIZE);
        ke_gas_density.data =
            0.5 * (state.gas.get_momentum_x().array().square() +
                   state.gas.get_momentum_y().array().square() +
                   state.gas.get_momentum_z().array().square());

        int total_cells =
            config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;
        for (int i = 0; i < total_cells; ++i) {
            if (state.gas.get_density().data[i] > 1e-12) {
                ke_gas_density.data[i] /= state.gas.get_density().data[i];
            } else {
                ke_gas_density.data[i] = 0.0;
            }
        }
        this->ke_gas = ke_gas_density.sum() * config.CELL_VOLUME;

        Grid3D internal_energy_density(config.MESH_SIZE);
        internal_energy_density.data =
            state.gas.get_pressure().array() / (config.GAMMA - 1.0);
        this->ie_gas = internal_energy_density.sum() * config.CELL_VOLUME;
    }
}