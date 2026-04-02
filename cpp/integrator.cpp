#include "integrator.h"
#include "particles.h"
#include "pocketfft_hdronly.h"
#include <chrono>

void update_cosmology( SimState& state, const Config& config ) {
    if( config.EXPANDING_UNIVERSE ) {
        double expansion_time = config.EXPANSION_START_T + state.total_time;
        state.scale_factor = std::pow( expansion_time, 2.0 / 3.0 );
        state.hubble_param = ( 2.0 / 3.0 ) / expansion_time;
    }
    else {
        state.scale_factor = 1.0;
        state.hubble_param = 0.0;
    }
}

Grid compute_gravitational_acceleration(
    GasGrid& gas,
    const Config& config,
    const Grid& dm_rho )
{
    Grid total_rho = dm_rho + ( config.USE_HYDRO ? gas.density : Grid::Zero( config.MESH_SIZE, config.MESH_SIZE ) );

    pocketfft::shape_t shape = { ( size_t )config.MESH_SIZE, ( size_t )config.MESH_SIZE };
    pocketfft::stride_t stride_r = { static_cast< ptrdiff_t >( sizeof( double ) * config.MESH_SIZE ), static_cast< ptrdiff_t >( sizeof( double ) ) };
    pocketfft::stride_t stride_c = { static_cast< ptrdiff_t >( sizeof( std::complex<double> ) * ( config.MESH_SIZE / 2 + 1 ) ), static_cast< ptrdiff_t >( sizeof( std::complex<double> ) ) };

    std::vector<std::complex<double>> rho_k( config.MESH_SIZE * ( config.MESH_SIZE / 2 + 1 ) );
    pocketfft::r2c( shape, stride_r, stride_c, { 0, 1 }, true, total_rho.data(), rho_k.data(), 1.0 );

    std::vector<std::complex<double>> phi_k( config.MESH_SIZE * ( config.MESH_SIZE / 2 + 1 ) );
    for( int i = 0; i < config.MESH_SIZE; ++i ) {
        for( int j = 0; j < config.MESH_SIZE / 2 + 1; ++j ) {
            if( i == 0 && j == 0 ) {
                phi_k[0] = { 0.0, 0.0 };
                continue;
            }
            double kx_freq = ( i < config.MESH_SIZE / 2 ) ? i : ( i - config.MESH_SIZE );
            double kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
            double ky = ( double )j * 2.0 * M_PI / config.DOMAIN_SIZE;
            double k2 = kx * kx + ky * ky;

            phi_k[i * ( config.MESH_SIZE / 2 + 1 ) + j] = rho_k[i * ( config.MESH_SIZE / 2 + 1 ) + j] * ( -2.0 * M_PI * config.G / sqrt( k2 ) );
        }
    }

    std::vector<double> phi_flat( config.MESH_SIZE * config.MESH_SIZE, 0.0 );
    pocketfft::c2r( shape, stride_c, stride_r, { 0, 1 }, false, phi_k.data(), phi_flat.data(), 1.0 );

    Grid phi = Eigen::Map<const Grid>( phi_flat.data(), config.MESH_SIZE, config.MESH_SIZE );
    double norm = 1.0 / ( config.MESH_SIZE * config.MESH_SIZE );

    gas.accel_x = (roll(phi, 1, 0) - roll(phi, -1, 0)) * norm / (2.0 * config.CELL_SIZE);
    gas.accel_y = (roll(phi, 1, 1) - roll(phi, -1, 1)) * norm / (2.0 * config.CELL_SIZE);

    return total_rho;
}

std::map<std::string, double> KDK_step(
    SimState& state,
    double dt,
    Config& config )
{
    std::map<std::string, double> timings;
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;

    update_cosmology( state, config );
    double a = state.scale_factor;
    double H = state.hubble_param;
    double a3 = a * a * a;

    // KICK 1 (t -> t + dt/2)
    state.gas.update_primitive_variables();
    Grid total_ax_gas = ( state.gas.accel_x.array() / a3 ) - ( 2 * H * state.gas.velocity_x.array() );
    Grid total_ay_gas = ( state.gas.accel_y.array() / a3 ) - ( 2 * H * state.gas.velocity_y.array() );
    Grid g_mom_x_source = state.gas.density.array() * total_ax_gas.array();
    Grid g_mom_y_source = state.gas.density.array() * total_ay_gas.array();
    Grid power_density = state.gas.velocity_x.array() * g_mom_x_source.array() + state.gas.velocity_y.array() * g_mom_y_source.array();

    if( config.USE_HYDRO ) {
        state.gas.momentum_x += g_mom_x_source * ( dt / 2.0 );
        state.gas.momentum_y += g_mom_y_source * ( dt / 2.0 );
        state.gas.energy += power_density * ( dt / 2.0 );
    }

    for( auto& p : state.dm.particles ) {
        double total_ax_p = ( p.acc.x / a3 ) - ( 2 * H * p.vel.x );
        double total_ay_p = ( p.acc.y / a3 ) - ( 2 * H * p.vel.y );
        p.vel.x += total_ax_p * dt / 2.0;
        p.vel.y += total_ay_p * dt / 2.0;
    }

    // DRIFT (t -> t + dt)
    for( auto& p : state.dm.particles ) {
        p.pos.x = fmod( p.pos.x + p.vel.x * dt + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
        p.pos.y = fmod( p.pos.y + p.vel.y * dt + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
    }

    start_time = std::chrono::high_resolution_clock::now();
    if( config.USE_HYDRO ) {
        state.gas.hydro_step( dt );
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["hydro"] = std::chrono::duration_cast< std::chrono::duration<double> >( end_time - start_time ).count();


    // UPDATE COSMOLOGY to t + dt
    state.total_time += dt;            // <-- Explicitly advance the clock
    update_cosmology( state, config ); // <-- Update the scale factor
    a = state.scale_factor;
    H = state.hubble_param;
    a3 = a * a * a;

    // COMPUTE FORCES at t + dt
    start_time = std::chrono::high_resolution_clock::now();

    // 1. Bin all particles ONCE
    state.dm.bin_and_assign_mass( config );

    // 2. Compute PM forces
    // (This function internally updates state.gas.accel_x and state.gas.accel_y)
    compute_gravitational_acceleration( state.gas, config, state.dm.dm_rho );

    std::vector<Vec2> pm_forces;
    if( config.USE_PM ) {
        state.dm.interpolate_cic_forces( state.gas.accel_x, state.gas.accel_y, pm_forces, config );
    }
    else {
        pm_forces.assign( state.dm.particles.size(), { 0.0, 0.0 } );
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["pm"] = std::chrono::duration_cast< std::chrono::duration<double> >( end_time - start_time ).count();


    start_time = std::chrono::high_resolution_clock::now();

    // 3. Compute PP forces
    std::vector<Vec2> pp_forces;
    if( config.USE_PP ) {
        state.dm.compute_pp_forces( pp_forces, config );
    }
    else {
        pp_forces.assign( state.dm.particles.size(), { 0.0, 0.0 } );
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["pp"] = std::chrono::duration_cast< std::chrono::duration<double> >( end_time - start_time ).count();

    // KICK 2 (t + dt/2 -> t + dt)
    state.gas.update_primitive_variables();
    total_ax_gas = ( state.gas.accel_x.array() / a3 ) - ( 2 * H * state.gas.velocity_x.array() );
    total_ay_gas = ( state.gas.accel_y.array() / a3 ) - ( 2 * H * state.gas.velocity_y.array() );
    g_mom_x_source = state.gas.density.array() * total_ax_gas.array();
    g_mom_y_source = state.gas.density.array() * total_ay_gas.array();
    power_density = state.gas.velocity_x.array() * g_mom_x_source.array() + state.gas.velocity_y.array() * g_mom_y_source.array();

    if( config.USE_HYDRO ) {
        state.gas.momentum_x += g_mom_x_source * ( dt / 2.0 );
        state.gas.momentum_y += g_mom_y_source * ( dt / 2.0 );
        state.gas.energy += power_density * ( dt / 2.0 );
    }

    for( size_t i = 0; i < state.dm.particles.size(); ++i ) {
        auto& p = state.dm.particles[i];
        Vec2 f = { pp_forces[i].x + pm_forces[i].x, pp_forces[i].y + pm_forces[i].y };
        p.acc.x = config.STANDING_PARTICLES ? 0 : f.x / p.mass;
        p.acc.y = config.STANDING_PARTICLES ? 0 : f.y / p.mass;
        double total_ax_p = ( p.acc.x / a3 ) - ( 2 * H * p.vel.x );
        double total_ay_p = ( p.acc.y / a3 ) - ( 2 * H * p.vel.y );
        p.vel.x += total_ax_p * dt / 2.0;
        p.vel.y += total_ay_p * dt / 2.0;
    }

    return timings;
}
