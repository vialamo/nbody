#include "ics.h"
#include "particles.h"
#include "integrator.h"
#include "pocketfft_hdronly.h"
#include <random>

static void create_zeldovich_ics( SimState& state, const Config& config ) {
    // Calculate the primordial fields at the Eulerian grid resolution (MESH_SIZE)
    int M = config.MESH_SIZE;
    double cell_size = config.DOMAIN_SIZE / M;

    size_t M3_real = static_cast< size_t >( M ) * M * M;
    size_t M3_complex = static_cast< size_t >( M ) * M * ( M / 2 + 1 );

    std::vector<double> real_space_random_field( M3_real );
    std::default_random_engine generator( config.SEED );
    std::normal_distribution<double> distribution( 0.0, 1.0 );
    for( auto& val : real_space_random_field ) {
        val = distribution( generator );
    }

    // FFT Setup
    pocketfft::shape_t shape_ic = { ( size_t )M, ( size_t )M, ( size_t )M };
    pocketfft::stride_t stride_r_ic = {
        static_cast< ptrdiff_t >( M * M * sizeof( double ) ),
        static_cast< ptrdiff_t >( M * sizeof( double ) ),
        sizeof( double )
    };
    pocketfft::stride_t stride_c_ic = {
        static_cast< ptrdiff_t >( M * ( M / 2 + 1 ) * sizeof( std::complex<double> ) ),
        static_cast< ptrdiff_t >( ( M / 2 + 1 ) * sizeof( std::complex<double> ) ),
        sizeof( std::complex<double> )
    };

    std::vector<std::complex<double>> random_k( M3_complex );
    pocketfft::r2c( shape_ic, stride_r_ic, stride_c_ic, { 0, 1, 2 }, true, real_space_random_field.data(), random_k.data(), 1.0 );

    std::vector<std::complex<double>> disp_x_k( M3_complex );
    std::vector<std::complex<double>> disp_y_k( M3_complex );
    std::vector<std::complex<double>> disp_z_k( M3_complex );

    for( int i = 0; i < M; ++i ) {
        for( int j = 0; j < M; ++j ) {
            for( int k = 0; k < M / 2 + 1; ++k ) {
                if( i == 0 && j == 0 && k == 0 ) {
                    disp_x_k[0] = { 0, 0 };
                    disp_y_k[0] = { 0, 0 };
                    disp_z_k[0] = { 0, 0 };
                    continue;
                }

                double kx_freq = static_cast< double >( ( i < M / 2 ) ? i : ( i - M ) );
                double ky_freq = static_cast< double >( ( j < M / 2 ) ? j : ( j - M ) );
                double kz_freq = static_cast< double >( k );

                double kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double ky = ky_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double kz = kz_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
                double k2 = kx * kx + ky * ky + kz * kz;

                size_t idx = static_cast< size_t >( i ) * M * ( M / 2 + 1 ) +
                    static_cast< size_t >( j ) * ( M / 2 + 1 ) +
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

    std::vector<double> disp_x_real( M3_real );
    std::vector<double> disp_y_real( M3_real );
    std::vector<double> disp_z_real( M3_real );

    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1, 2 }, false, disp_x_k.data(), disp_x_real.data(), 1.0 );
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1, 2 }, false, disp_y_k.data(), disp_y_real.data(), 1.0 );
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1, 2 }, false, disp_z_k.data(), disp_z_real.data(), 1.0 );

    double norm_ic = 1.0 / static_cast< double >( M3_real );
    double std_x = 0, std_y = 0, std_z = 0;
    for( size_t i = 0; i < M3_real; ++i ) {
        disp_x_real[i] *= norm_ic;
        disp_y_real[i] *= norm_ic;
        disp_z_real[i] *= norm_ic;
        std_x += disp_x_real[i] * disp_x_real[i];
        std_y += disp_y_real[i] * disp_y_real[i];
        std_z += disp_z_real[i] * disp_z_real[i];
    }
    std_x = sqrt( std_x / M3_real );
    std_y = sqrt( std_y / M3_real );
    std_z = sqrt( std_z / M3_real );

    // Populate Gas Velocities directly onto the MESH_SIZE grid
    if( config.USE_HYDRO ) {
        for( size_t i = 0; i < M3_real; ++i ) {
            double dx = ( disp_x_real[i] / std_x ) * state.scale_factor * cell_size;
            double dy = ( disp_y_real[i] / std_y ) * state.scale_factor * cell_size;
            double dz = ( disp_z_real[i] / std_z ) * state.scale_factor * cell_size;

            state.gas.velocity_x.data[i] = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dx;
            state.gas.velocity_y.data[i] = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dy;
            state.gas.velocity_z.data[i] = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dz;
        }
    }

    // Map Particles from the High-Res Grid
    state.dm.particles.clear();
    int N_part = config.N_PER_SIDE;
    double spacing = config.DOMAIN_SIZE / N_part;

    for( int i = 0; i < N_part; ++i ) {
        for( int j = 0; j < N_part; ++j ) {
            for( int k = 0; k < N_part; ++k ) {
                // Particle's perfectly uniform starting coordinate
                double qx = ( i + 0.5 ) * spacing;
                double qy = ( j + 0.5 ) * spacing;
                double qz = ( k + 0.5 ) * spacing;

                // Map physical coordinate to the high-res MESH_SIZE grid index
                int ix = static_cast< int >( qx / cell_size ) % M;
                int iy = static_cast< int >( qy / cell_size ) % M;
                int iz = static_cast< int >( qz / cell_size ) % M;

                size_t idx = static_cast< size_t >( ix ) * M * M +
                    static_cast< size_t >( iy ) * M +
                    static_cast< size_t >( iz );

                double dx = ( disp_x_real[idx] / std_x ) * state.scale_factor * cell_size;
                double dy = ( disp_y_real[idx] / std_y ) * state.scale_factor * cell_size;
                double dz = ( disp_z_real[idx] / std_z ) * state.scale_factor * cell_size;

                Particle p;
                p.pos.x = fmod( qx + dx + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
                p.pos.y = fmod( qy + dy + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
                p.pos.z = fmod( qz + dz + config.DOMAIN_SIZE, config.DOMAIN_SIZE );

                p.vel.x = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dx;
                p.vel.y = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dy;
                p.vel.z = config.STANDING_PARTICLES ? 0.0 : state.hubble_param * dz;
                p.mass = config.DM_PARTICLE_MASS;

                state.dm.particles.push_back( p );
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

    state.dm.bin_and_assign_mass( config );

    // Initialize Gas to perfectly trace the Dark Matter
    if( config.USE_HYDRO ) {
        double total_dm_mass = state.dm.particles.size() * config.DM_PARTICLE_MASS;
        double mass_ratio = config.GAS_TOTAL_MASS / total_dm_mass;

        state.gas.density.data = state.dm.dm_rho.data * mass_ratio;

        // Prevent hydro crashes in deep voids by setting a density floor
        state.gas.density.data = ( state.gas.density.array() < 1e-12 ).select( 1e-12, state.gas.density.data );

        state.gas.momentum_x.data = state.gas.density.array() * state.gas.velocity_x.array();
        state.gas.momentum_y.data = state.gas.density.array() * state.gas.velocity_y.array();
        state.gas.momentum_z.data = state.gas.density.array() * state.gas.velocity_z.array();

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