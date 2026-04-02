#include "ics.h"
#include "particles.h"
#include "integrator.h"
#include "pocketfft_hdronly.h"
#include <random>

static void create_zeldovich_ics( SimState& state, const Config& config ) {
    int n_per_side = config.N_PER_SIDE;
    int ic_mesh_size = n_per_side;
    double cell_size = config.DOMAIN_SIZE / ic_mesh_size;

    std::vector<double> real_space_random_field( ic_mesh_size * ic_mesh_size );
    std::default_random_engine generator( config.SEED );
    std::normal_distribution<double> distribution( 0.0, 1.0 );
    for( auto& val : real_space_random_field ) {
        val = distribution( generator );
    }

    pocketfft::shape_t shape_ic = { ( size_t )ic_mesh_size, ( size_t )ic_mesh_size };
    pocketfft::stride_t stride_r_ic = { static_cast< ptrdiff_t >( sizeof( double ) * ic_mesh_size ), sizeof( double ) };
    pocketfft::stride_t stride_c_ic = { static_cast< ptrdiff_t >( sizeof( std::complex<double> ) ) * ( ic_mesh_size / 2 + 1 ), sizeof( std::complex<double> ) };

    std::vector<std::complex<double>> random_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );
    pocketfft::r2c( shape_ic, stride_r_ic, stride_c_ic, { 0, 1 }, true, real_space_random_field.data(), random_k.data(), 1.0 );

    std::vector<std::complex<double>> disp_x_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );
    std::vector<std::complex<double>> disp_y_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );

    for( int i = 0; i < ic_mesh_size; ++i ) {
        for( int j = 0; j < ic_mesh_size / 2 + 1; ++j ) {
            double kx_freq = ( i < ic_mesh_size / 2 ) ? i : ( i - ic_mesh_size );
            double ky_freq = j;
            double kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
            double ky = ky_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
            double k2 = kx * kx + ky * ky;
            int idx = i * ( ic_mesh_size / 2 + 1 ) + j;
            if( i == 0 && j == 0 ) {
                disp_x_k[idx] = { 0, 0 };
                disp_y_k[idx] = { 0, 0 };
                continue;
            }
            double power_spectrum_sqrt = sqrt( pow( k2, config.INITIAL_POWER_SPECTRUM_INDEX / 2.0 ) );
            std::complex<double> delta_k = random_k[idx] * power_spectrum_sqrt;
            std::complex<double> phi_k = -delta_k / k2;
            disp_x_k[idx] = std::complex<double>( 0, -1 ) * kx * phi_k;
            disp_y_k[idx] = std::complex<double>( 0, -1 ) * ky * phi_k;
        }
    }

    std::vector<double> disp_x_real( ic_mesh_size * ic_mesh_size );
    std::vector<double> disp_y_real( ic_mesh_size * ic_mesh_size );
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1 }, false, disp_x_k.data(), disp_x_real.data(), 1.0 );
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1 }, false, disp_y_k.data(), disp_y_real.data(), 1.0 );

    double norm_ic = 1.0 / ( ic_mesh_size * ic_mesh_size );
    double std_x = 0, std_y = 0;
    for( size_t i = 0; i < disp_x_real.size(); ++i ) {
        disp_x_real[i] *= norm_ic;
        disp_y_real[i] *= norm_ic;
        std_x += disp_x_real[i] * disp_x_real[i];
        std_y += disp_y_real[i] * disp_y_real[i];
    }
    std_x = sqrt( std_x / disp_x_real.size() );
    std_y = sqrt( std_y / disp_y_real.size() );

    state.dm.particles.clear();
    double spacing = config.DOMAIN_SIZE / n_per_side;
    for( int i = 0; i < n_per_side; ++i ) {
        for( int j = 0; j < n_per_side; ++j ) {
            double qx = ( i + 0.5 ) * spacing;
            double qy = ( j + 0.5 ) * spacing;
            double dx = ( disp_x_real[i * n_per_side + j] / std_x ) * state.scale_factor * cell_size;
            double dy = ( disp_y_real[i * n_per_side + j] / std_y ) * state.scale_factor * cell_size;
            Particle p;
            p.pos.x = fmod( qx + dx + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
            p.pos.y = fmod( qy + dy + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
            p.vel.x = config.STANDING_PARTICLES ? 0 : state.hubble_param * dx;
            p.vel.y = config.STANDING_PARTICLES ? 0 : state.hubble_param * dy;
            p.mass = config.DM_PARTICLE_MASS;
            state.dm.particles.push_back( p );
        }
    }
}

SimState initialize_state( Config& config ) {
    // Create the empty state object
    SimState state( config );

    // Populate Particles
    state.total_time = 0;
    update_cosmology( state, config );

    create_zeldovich_ics( state, config );

    // Populate Gas
    if( config.USE_HYDRO ) {
        state.gas.density.fill( config.GAS_TOTAL_MASS / ( config.DOMAIN_SIZE * config.DOMAIN_SIZE ) );
        double initial_internal_energy = 1e-6;
        state.gas.energy = state.gas.density.array() * initial_internal_energy;
        state.gas.update_primitive_variables();
    }

    // Compute Step 0 Forces
    state.dm.bin_and_assign_mass( config );
    Grid total_rho = compute_gravitational_acceleration( state.gas, config, state.dm.dm_rho );

    std::vector<Vec2> pp_forces, pm_forces;
    if( config.USE_PM ) { state.dm.interpolate_cic_forces( state.gas.accel_x, state.gas.accel_y, pm_forces, config ); }
    else { pm_forces.assign( state.dm.particles.size(), { 0.0, 0.0 } ); }

    if( config.USE_PP ) { state.dm.compute_pp_forces( pp_forces, config ); }
    else { pp_forces.assign( state.dm.particles.size(), { 0.0, 0.0 } ); }

    for( size_t i = 0; i < state.dm.particles.size(); ++i ) {
        state.dm.particles[i].acc.x = ( pp_forces[i].x + pm_forces[i].x ) / state.dm.particles[i].mass;
        state.dm.particles[i].acc.y = ( pp_forces[i].y + pm_forces[i].y ) / state.dm.particles[i].mass;
    }

    return state;
}