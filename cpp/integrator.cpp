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

Grid3D compute_gravitational_acceleration( GasGrid& gas, const Config& config, const Grid3D& dm_rho ) {
    int N = config.MESH_SIZE;
    Grid3D total_rho( N );
    total_rho.data = dm_rho.data + ( config.USE_HYDRO ? gas.density.data : Eigen::VectorXd::Zero( N * N * N ) );

    pocketfft::shape_t shape = { ( size_t )N, ( size_t )N, ( size_t )N };

    pocketfft::stride_t stride_r = {
        static_cast< ptrdiff_t >( ( size_t )N * N * sizeof( double ) ),
        static_cast< ptrdiff_t >( N * sizeof( double ) ),
        sizeof( double )
    };
    pocketfft::stride_t stride_c = {
        static_cast< ptrdiff_t >( ( size_t )N * ( N / 2 + 1 ) * sizeof( std::complex<double> ) ),
        static_cast< ptrdiff_t >( ( N / 2 + 1 ) * sizeof( std::complex<double> ) ),
        sizeof( std::complex<double> )
    };

    std::vector<std::complex<double>> rho_k( ( size_t )N * N * ( N / 2 + 1 ) );
    pocketfft::r2c( shape, stride_r, stride_c, { 0, 1, 2 }, true, total_rho.raw_data(), rho_k.data(), 1.0 );

    std::vector<std::complex<double>> phi_k( ( size_t )N * N * ( N / 2 + 1 ) );

    for( int i = 0; i < N; ++i ) {
        for( int j = 0; j < N; ++j ) {
            for( int k = 0; k < N / 2 + 1; ++k ) {
                if( i == 0 && j == 0 && k == 0 ) {
                    phi_k[0] = { 0.0, 0.0 };
                    continue;
                }

                double kx_freq = static_cast< double >( ( i < N / 2 ) ? i : ( i - N ) );
                double ky_freq = static_cast< double >( ( j < N / 2 ) ? j : ( j - N ) );
                double kz_freq = static_cast< double >( k );

                double kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double ky = ky_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double kz = kz_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double k2 = kx * kx + ky * ky + kz * kz;

                size_t idx = static_cast< size_t >( i ) * N * ( N / 2 + 1 ) +
                    static_cast< size_t >( j ) * ( N / 2 + 1 ) +
                    static_cast< size_t >( k );

                phi_k[idx] = rho_k[idx] * ( -4.0 * M_PI * config.G / k2 );
            }
        }
    }

    Grid3D phi( N );
    pocketfft::c2r( shape, stride_c, stride_r, { 0, 1, 2 }, false, phi_k.data(), const_cast< double* >( phi.raw_data() ), 1.0 );

    double norm = 1.0 / ( ( double )N * N * N );
    double factor = -norm / ( 2.0 * config.CELL_SIZE );

    for( int i = 0; i < N; ++i ) {
        for( int j = 0; j < N; ++j ) {
            for( int k = 0; k < N; ++k ) {
                gas.accel_x( i, j, k ) = ( phi( ( i + 1 ) % N, j, k ) - phi( ( i - 1 + N ) % N, j, k ) ) * factor;
                gas.accel_y( i, j, k ) = ( phi( i, ( j + 1 ) % N, k ) - phi( i, ( j - 1 + N ) % N, k ) ) * factor;
                gas.accel_z( i, j, k ) = ( phi( i, j, ( k + 1 ) % N ) - phi( i, j, ( k - 1 + N ) % N ) ) * factor;
            }
        }
    }

    return total_rho;
}

std::map<std::string, double> KDK_step( SimState& state, double dt, Config& config ) {
    std::map<std::string, double> timings;
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;

    update_cosmology( state, config );
    double a = state.scale_factor;
    double H = state.hubble_param;
    double a3 = a * a * a;

    // KICK 1
    state.gas.update_primitive_variables();
    Grid3D total_ax_gas( config.MESH_SIZE ), total_ay_gas( config.MESH_SIZE ), total_az_gas( config.MESH_SIZE );
    total_ax_gas.data = ( state.gas.accel_x.array() / a3 ) - ( 2 * H * state.gas.velocity_x.array() );
    total_ay_gas.data = ( state.gas.accel_y.array() / a3 ) - ( 2 * H * state.gas.velocity_y.array() );
    total_az_gas.data = ( state.gas.accel_z.array() / a3 ) - ( 2 * H * state.gas.velocity_z.array() );

    Grid3D g_mom_x_source( config.MESH_SIZE ), g_mom_y_source( config.MESH_SIZE ), g_mom_z_source( config.MESH_SIZE );
    g_mom_x_source.data = state.gas.density.array() * total_ax_gas.array();
    g_mom_y_source.data = state.gas.density.array() * total_ay_gas.array();
    g_mom_z_source.data = state.gas.density.array() * total_az_gas.array();

    Grid3D power_density( config.MESH_SIZE );
    power_density.data = state.gas.velocity_x.array() * g_mom_x_source.array() +
        state.gas.velocity_y.array() * g_mom_y_source.array() +
        state.gas.velocity_z.array() * g_mom_z_source.array();

    if( config.USE_HYDRO ) {
        state.gas.momentum_x.array() += g_mom_x_source.array() * ( dt / 2.0 );
        state.gas.momentum_y.array() += g_mom_y_source.array() * ( dt / 2.0 );
        state.gas.momentum_z.array() += g_mom_z_source.array() * ( dt / 2.0 );
        state.gas.energy.array() += power_density.array() * ( dt / 2.0 );
    }

    for( auto& p : state.dm.particles ) {
        double total_ax_p = ( p.acc.x / a3 ) - ( 2 * H * p.vel.x );
        double total_ay_p = ( p.acc.y / a3 ) - ( 2 * H * p.vel.y );
        double total_az_p = ( p.acc.z / a3 ) - ( 2 * H * p.vel.z );
        p.vel.x += total_ax_p * dt / 2.0;
        p.vel.y += total_ay_p * dt / 2.0;
        p.vel.z += total_az_p * dt / 2.0;
    }

    // DRIFT
    for( auto& p : state.dm.particles ) {
        p.pos.x = fmod( p.pos.x + p.vel.x * dt + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
        p.pos.y = fmod( p.pos.y + p.vel.y * dt + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
        p.pos.z = fmod( p.pos.z + p.vel.z * dt + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
    }

    start_time = std::chrono::high_resolution_clock::now();
    if( config.USE_HYDRO ) {
        state.gas.hydro_step( dt );
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["hydro"] = std::chrono::duration_cast< std::chrono::duration<double> >( end_time - start_time ).count();

    // UPDATE COSMOLOGY to t + dt
    state.total_time += dt;
    update_cosmology( state, config );
    a = state.scale_factor;
    H = state.hubble_param;
    a3 = a * a * a;

    // COMPUTE FORCES
    start_time = std::chrono::high_resolution_clock::now();
    state.dm.bin_and_assign_mass( config );
    compute_gravitational_acceleration( state.gas, config, state.dm.dm_rho );

    std::vector<Vec3> pm_forces;
    if( config.USE_PM ) {
        state.dm.interpolate_cic_forces( state.gas.accel_x, state.gas.accel_y, state.gas.accel_z, pm_forces, config );
    }
    else {
        pm_forces.assign( state.dm.particles.size(), { 0.0, 0.0, 0.0 } );
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["pm"] = std::chrono::duration_cast< std::chrono::duration<double> >( end_time - start_time ).count();

    start_time = std::chrono::high_resolution_clock::now();
    std::vector<Vec3> pp_forces;
    if( config.USE_PP ) {
        state.dm.compute_pp_forces( pp_forces, config );
    }
    else {
        pp_forces.assign( state.dm.particles.size(), { 0.0, 0.0, 0.0 } );
    }
    end_time = std::chrono::high_resolution_clock::now();
    timings["pp"] = std::chrono::duration_cast< std::chrono::duration<double> >( end_time - start_time ).count();

    // KICK 2
    state.gas.update_primitive_variables();
    total_ax_gas.data = ( state.gas.accel_x.array() / a3 ) - ( 2 * H * state.gas.velocity_x.array() );
    total_ay_gas.data = ( state.gas.accel_y.array() / a3 ) - ( 2 * H * state.gas.velocity_y.array() );
    total_az_gas.data = ( state.gas.accel_z.array() / a3 ) - ( 2 * H * state.gas.velocity_z.array() );

    g_mom_x_source.data = state.gas.density.array() * total_ax_gas.array();
    g_mom_y_source.data = state.gas.density.array() * total_ay_gas.array();
    g_mom_z_source.data = state.gas.density.array() * total_az_gas.array();

    power_density.data = state.gas.velocity_x.array() * g_mom_x_source.array() +
        state.gas.velocity_y.array() * g_mom_y_source.array() +
        state.gas.velocity_z.array() * g_mom_z_source.array();

    if( config.USE_HYDRO ) {
        state.gas.momentum_x.array() += g_mom_x_source.array() * ( dt / 2.0 );
        state.gas.momentum_y.array() += g_mom_y_source.array() * ( dt / 2.0 );
        state.gas.momentum_z.array() += g_mom_z_source.array() * ( dt / 2.0 );
        state.gas.energy.array() += power_density.array() * ( dt / 2.0 );
    }

    for( size_t i = 0; i < state.dm.particles.size(); ++i ) {
        auto& p = state.dm.particles[i];
        Vec3 f = { pp_forces[i].x + pm_forces[i].x, pp_forces[i].y + pm_forces[i].y, pp_forces[i].z + pm_forces[i].z };
        p.acc.x = config.STANDING_PARTICLES ? 0 : f.x / p.mass;
        p.acc.y = config.STANDING_PARTICLES ? 0 : f.y / p.mass;
        p.acc.z = config.STANDING_PARTICLES ? 0 : f.z / p.mass;
        p.vel.x += ( ( p.acc.x / a3 ) - ( 2 * H * p.vel.x ) ) * dt / 2.0;
        p.vel.y += ( ( p.acc.y / a3 ) - ( 2 * H * p.vel.y ) ) * dt / 2.0;
        p.vel.z += ( ( p.acc.z / a3 ) - ( 2 * H * p.vel.z ) ) * dt / 2.0;
    }

    return timings;
}