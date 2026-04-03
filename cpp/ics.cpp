#include "ics.h"
#include "particles.h"
#include "integrator.h"
#include "pocketfft_hdronly.h"
#include <random>

static void create_zeldovich_ics( SimState& state, const Config& config ) {
    int N = config.N_PER_SIDE;
    double cell_size = config.DOMAIN_SIZE / N;

    size_t N3_real = static_cast< size_t >( N ) * N * N;
    size_t N3_complex = static_cast< size_t >( N ) * N * ( N / 2 + 1 );

    std::vector<double> real_space_random_field( N3_real );
    std::default_random_engine generator( config.SEED );
    std::normal_distribution<double> distribution( 0.0, 1.0 );
    for( auto& val : real_space_random_field ) {
        val = distribution( generator );
    }

    // FFT Setup
    pocketfft::shape_t shape_ic = { ( size_t )N, ( size_t )N, ( size_t )N };
    pocketfft::stride_t stride_r_ic = {
        static_cast< ptrdiff_t >( N * N * sizeof( double ) ),
        static_cast< ptrdiff_t >( N * sizeof( double ) ),
        sizeof( double )
    };
    pocketfft::stride_t stride_c_ic = {
        static_cast< ptrdiff_t >( N * ( N / 2 + 1 ) * sizeof( std::complex<double> ) ),
        static_cast< ptrdiff_t >( ( N / 2 + 1 ) * sizeof( std::complex<double> ) ),
        sizeof( std::complex<double> )
    };

    std::vector<std::complex<double>> random_k( N3_complex );
    pocketfft::r2c( shape_ic, stride_r_ic, stride_c_ic, { 0, 1, 2 }, true, real_space_random_field.data(), random_k.data(), 1.0 );

    std::vector<std::complex<double>> disp_x_k( N3_complex );
    std::vector<std::complex<double>> disp_y_k( N3_complex );
    std::vector<std::complex<double>> disp_z_k( N3_complex );

    for( int i = 0; i < N; ++i ) {
        for( int j = 0; j < N; ++j ) {
            for( int k = 0; k < N / 2 + 1; ++k ) {
                if( i == 0 && j == 0 && k == 0 ) {
                    disp_x_k[0] = { 0, 0 };
                    disp_y_k[0] = { 0, 0 };
                    disp_z_k[0] = { 0, 0 };
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

                double power_spectrum_sqrt = sqrt( pow( k2, config.INITIAL_POWER_SPECTRUM_INDEX / 2.0 ) );
                std::complex<double> delta_k = random_k[idx] * power_spectrum_sqrt;
                std::complex<double> phi_k = -delta_k / k2;

                disp_x_k[idx] = std::complex<double>( 0, -1 ) * kx * phi_k;
                disp_y_k[idx] = std::complex<double>( 0, -1 ) * ky * phi_k;
                disp_z_k[idx] = std::complex<double>( 0, -1 ) * kz * phi_k;
            }
        }
    }

    std::vector<double> disp_x_real( N3_real );
    std::vector<double> disp_y_real( N3_real );
    std::vector<double> disp_z_real( N3_real );

    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1, 2 }, false, disp_x_k.data(), disp_x_real.data(), 1.0 );
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1, 2 }, false, disp_y_k.data(), disp_y_real.data(), 1.0 );
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1, 2 }, false, disp_z_k.data(), disp_z_real.data(), 1.0 );

    double norm_ic = 1.0 / static_cast< double >( N3_real );
    double std_x = 0, std_y = 0, std_z = 0;
    for( size_t i = 0; i < N3_real; ++i ) {
        disp_x_real[i] *= norm_ic;
        disp_y_real[i] *= norm_ic;
        disp_z_real[i] *= norm_ic;
        std_x += disp_x_real[i] * disp_x_real[i];
        std_y += disp_y_real[i] * disp_y_real[i];
        std_z += disp_z_real[i] * disp_z_real[i];
    }
    std_x = sqrt( std_x / N3_real );
    std_y = sqrt( std_y / N3_real );
    std_z = sqrt( std_z / N3_real );

    state.dm.particles.clear();
    double spacing = config.DOMAIN_SIZE / N;

    for( int i = 0; i < N; ++i ) {
        for( int j = 0; j < N; ++j ) {
            for( int k = 0; k < N; ++k ) {
                double qx = ( i + 0.5 ) * spacing;
                double qy = ( j + 0.5 ) * spacing;
                double qz = ( k + 0.5 ) * spacing;

                size_t idx = static_cast< size_t >( i ) * N * N + static_cast< size_t >( j ) * N + static_cast< size_t >( k );

                double dx = ( disp_x_real[idx] / std_x ) * state.scale_factor * cell_size;
                double dy = ( disp_y_real[idx] / std_y ) * state.scale_factor * cell_size;
                double dz = ( disp_z_real[idx] / std_z ) * state.scale_factor * cell_size;

                Particle p;
                p.pos.x = fmod( qx + dx + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
                p.pos.y = fmod( qy + dy + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
                p.pos.z = fmod( qz + dz + config.DOMAIN_SIZE, config.DOMAIN_SIZE );

                p.vel.x = config.STANDING_PARTICLES ? 0 : state.hubble_param * dx;
                p.vel.y = config.STANDING_PARTICLES ? 0 : state.hubble_param * dy;
                p.vel.z = config.STANDING_PARTICLES ? 0 : state.hubble_param * dz;
                p.mass = config.DM_PARTICLE_MASS;

                state.dm.particles.push_back( p );

                // Save the primordial velocities to the Gas grid
                if( config.USE_HYDRO ) {
                    size_t flat_idx = static_cast< size_t >( i ) * N * N + static_cast< size_t >( j ) * N + static_cast< size_t >( k );
                    state.gas.velocity_x.data[flat_idx] = p.vel.x;
                    state.gas.velocity_y.data[flat_idx] = p.vel.y;
                    state.gas.velocity_z.data[flat_idx] = p.vel.z;
                }
            }
        }
    }
}

SimState initialize_state( Config& config ) {
    SimState state( config );
    state.total_time = 0;
    update_cosmology( state, config );

    // Generate Particle Positions and Velocities
    create_zeldovich_ics( state, config );

    // Bin Dark Matter to get the primordial density structure
    state.dm.bin_and_assign_mass( config );

    // Initialize Gas to perfectly trace the Dark Matter
    if( config.USE_HYDRO ) {
        // Gas density is directly proportional to DM density
        double total_dm_mass = state.dm.particles.size() * config.DM_PARTICLE_MASS;
        double mass_ratio = config.GAS_TOTAL_MASS / total_dm_mass;

        state.gas.density.data = state.dm.dm_rho.data * mass_ratio;

        // Prevent hydro crashes in deep voids by setting a density floor
        state.gas.density.data = ( state.gas.density.array() < 1e-12 ).select( 1e-12, state.gas.density.data );

        // Momentum = Density * Velocity
        state.gas.momentum_x.data = state.gas.density.array() * state.gas.velocity_x.array();
        state.gas.momentum_y.data = state.gas.density.array() * state.gas.velocity_y.array();
        state.gas.momentum_z.data = state.gas.density.array() * state.gas.velocity_z.array();

        // Total Energy = Internal + Kinetic
        double initial_internal_energy = 1e-6;
        Grid3D kin_energy( config.MESH_SIZE );
        kin_energy.data = 0.5 * ( state.gas.momentum_x.array().square() +
            state.gas.momentum_y.array().square() +
            state.gas.momentum_z.array().square() ) / state.gas.density.array();

        state.gas.energy.data = ( state.gas.density.array() * initial_internal_energy ) + kin_energy.data.array();
        state.gas.update_primitive_variables();
    }

    // Compute Step 0 Gravity Forces
    Grid3D total_rho = compute_gravitational_acceleration( state.gas, config, state.dm.dm_rho );

    std::vector<Vec3> pp_forces, pm_forces;
    if( config.USE_PM ) { state.dm.interpolate_cic_forces( state.gas.accel_x, state.gas.accel_y, state.gas.accel_z, pm_forces, config ); }
    else { pm_forces.assign( state.dm.particles.size(), { 0.0, 0.0, 0.0 } ); }

    if( config.USE_PP ) { state.dm.compute_pp_forces( pp_forces, config ); }
    else { pp_forces.assign( state.dm.particles.size(), { 0.0, 0.0, 0.0 } ); }

    for( size_t i = 0; i < state.dm.particles.size(); ++i ) {
        state.dm.particles[i].acc.x = ( pp_forces[i].x + pm_forces[i].x ) / state.dm.particles[i].mass;
        state.dm.particles[i].acc.y = ( pp_forces[i].y + pm_forces[i].y ) / state.dm.particles[i].mass;
        state.dm.particles[i].acc.z = ( pp_forces[i].z + pm_forces[i].z ) / state.dm.particles[i].mass;
    }

    return state;
}