#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <complex>
#include <chrono>

// Include the SFML graphics header
#include <SFML/Graphics.hpp>

// Include the header-only pocketfft library
#include "pocketfft_hdronly.h"

// Define M_PI if it's not provided by the compiler (it's non-standard)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ------------------------
// Vector and Particle Structs
// ------------------------
struct Vec2 {
    double x = 0.0, y = 0.0;
};

struct Particle {
    Vec2 pos;
    Vec2 vel;
    Vec2 acc;
    double mass;
};

// ------------------------
// Simulation parameters
// ------------------------
// Domain and mesh
const double DOMAIN_SIZE = 1.0;
const int MESH_SIZE = 16;
const double CELL_SIZE = DOMAIN_SIZE / MESH_SIZE;

// P3M parameters
const bool USE_PM = true;
const double CUTOFF_RADIUS = 2.5 * CELL_SIZE;
const double CUTOFF_RADIUS_SQUARED = CUTOFF_RADIUS * CUTOFF_RADIUS;
const double CUTOFF_TRANSITION_WIDTH = 0.2 * CUTOFF_RADIUS;
const double R_SWITCH_START = CUTOFF_RADIUS - CUTOFF_TRANSITION_WIDTH;
const double R_SWITCH_START_SQ = R_SWITCH_START * R_SWITCH_START;

// Particles
const int N_PER_SIDE = 50;
const int NUM_PARTICLES = N_PER_SIDE * N_PER_SIDE;
const double TOTAL_MASS = 1.0;
const double PARTICLE_MASS = 1.0 / NUM_PARTICLES;

// Physics
const double MEAN_INTERPARTICLE_SPACING = DOMAIN_SIZE / sqrt( NUM_PARTICLES );
const double SOFTENING_SQUARED = pow( MEAN_INTERPARTICLE_SPACING / 50.0, 2 );
const bool EXPANDING_UNIVERSE = true;
const double EXPANSION_START_T = 0.1;
const double G = 1.0 / ( 6.0 * M_PI );

// Time integration
const double DYNAMICAL_TIME = 1.0 / sqrt( G );
const double DT = 1e-5 * DYNAMICAL_TIME;


// Misc
const int DEBUG_INFO_EVERY_CYCLES = 40;
const int SAVE_RENDER_EVERY_CYCLES = 0;
const int SEED = 42;

// Render
const unsigned int RENDER_SIZE = 800;
const double RENDER_SCALE = RENDER_SIZE / DOMAIN_SIZE;


// Global cosmological variables
double scale_factor = 1.0;
double hubble_param = 0.0;


// ------------------------
// Helper Functions
// ------------------------
double displacement( double dx ) {
    return fmod( dx + 0.5 * DOMAIN_SIZE + DOMAIN_SIZE, DOMAIN_SIZE ) - 0.5 * DOMAIN_SIZE;
}

void update_cosmology( double total_time ) {
    if( EXPANDING_UNIVERSE ) {
        double expansion_time = EXPANSION_START_T + total_time;
        scale_factor = pow( expansion_time, 2.0 / 3.0 );
        hubble_param = ( 2.0 / 3.0 ) / expansion_time;
    }
    else {
        scale_factor = 1.0;
        hubble_param = 0.0;
    }
}

// ------------------------
// Force Calculation
// ------------------------
struct CIC_Data {
    int ix, iy;
    double w1, w2, w3, w4;
};

void compute_mesh_forces( std::vector<Particle>& particles, std::vector<Vec2>& mesh_forces ) {
    std::vector<double> rho( MESH_SIZE * MESH_SIZE, 0.0 );
    std::vector<CIC_Data> cic_data( particles.size() );

    // 1. Mass Assignment (CIC)
    for( size_t i = 0; i < particles.size(); ++i ) {
        const auto& p = particles[i];
        int ix = static_cast< int >( p.pos.x / CELL_SIZE );
        int iy = static_cast< int >( p.pos.y / CELL_SIZE );
        double frac_x = ( p.pos.x / CELL_SIZE ) - ix;
        double frac_y = ( p.pos.y / CELL_SIZE ) - iy;

        double w1 = ( 1 - frac_x ) * ( 1 - frac_y );
        double w2 = frac_x * ( 1 - frac_y );
        double w3 = ( 1 - frac_x ) * frac_y;
        double w4 = frac_x * frac_y;

        cic_data[i] = { ix, iy, w1, w2, w3, w4 };

        rho[( ix % MESH_SIZE ) * MESH_SIZE + ( iy % MESH_SIZE )] += p.mass * w1;
        rho[( ( ix + 1 ) % MESH_SIZE ) * MESH_SIZE + ( iy % MESH_SIZE )] += p.mass * w2;
        rho[( ix % MESH_SIZE ) * MESH_SIZE + ( ( iy + 1 ) % MESH_SIZE )] += p.mass * w3;
        rho[( ( ix + 1 ) % MESH_SIZE ) * MESH_SIZE + ( ( iy + 1 ) % MESH_SIZE )] += p.mass * w4;
    }

    for( auto& val : rho ) {
        val /= ( CELL_SIZE * CELL_SIZE );
    }

    // 2. FFT and Poisson Solver using pocketfft
    pocketfft::shape_t shape = { ( size_t )MESH_SIZE, ( size_t )MESH_SIZE };
    pocketfft::stride_t stride_r = { sizeof( double ) * MESH_SIZE, sizeof( double ) };
    pocketfft::stride_t stride_c = { sizeof( std::complex<double> ) * ( MESH_SIZE / 2 + 1 ), sizeof( std::complex<double> ) };

    std::vector<std::complex<double>> rho_k( MESH_SIZE * ( MESH_SIZE / 2 + 1 ) );
    pocketfft::r2c( shape, stride_r, stride_c, { 0, 1 }, true, rho.data(), rho_k.data(), 1.0 );


    std::vector<std::complex<double>> phi_k( MESH_SIZE * ( MESH_SIZE / 2 + 1 ) );
    for( int i = 0; i < MESH_SIZE; ++i ) {
        for( int j = 0; j < MESH_SIZE / 2 + 1; ++j ) {
            double kx_freq = ( i < MESH_SIZE / 2 ) ? i : ( i - MESH_SIZE );
            double ky_freq = j;

            double kx = kx_freq * 2.0 * M_PI / DOMAIN_SIZE;
            double ky = ky_freq * 2.0 * M_PI / DOMAIN_SIZE;
            double k2 = kx * kx + ky * ky;
            double k = sqrt( k2 );

            int idx = i * ( MESH_SIZE / 2 + 1 ) + j;
            if( i == 0 && j == 0 ) {
                phi_k[idx] = { 0.0, 0.0 };
                continue;
            }

            double factor = -2.0 * M_PI * G / k;
            phi_k[idx] = rho_k[idx] * factor;
        }
    }

    // 3. Inverse FFT
    std::vector<double> phi( MESH_SIZE * MESH_SIZE, 0.0 );
    pocketfft::c2r( shape, stride_c, stride_r, { 0, 1 }, false, phi_k.data(), phi.data(), 1.0 );

    // 4. Gradient
    std::vector<Vec2> force_grid( MESH_SIZE * MESH_SIZE );
    double norm = 1.0 / ( MESH_SIZE * MESH_SIZE ); // Normalization for pocketfft
    for( int i = 0; i < MESH_SIZE; ++i ) {
        for( int j = 0; j < MESH_SIZE; ++j ) {
            int current = i * MESH_SIZE + j;
            int prev_x = ( ( i - 1 + MESH_SIZE ) % MESH_SIZE ) * MESH_SIZE + j;
            int next_x = ( ( i + 1 ) % MESH_SIZE ) * MESH_SIZE + j;
            int prev_y = i * MESH_SIZE + ( ( j - 1 + MESH_SIZE ) % MESH_SIZE );
            int next_y = i * MESH_SIZE + ( ( j + 1 ) % MESH_SIZE );

            force_grid[current].x = ( phi[prev_x] - phi[next_x] ) * norm / ( 2.0 * CELL_SIZE );
            force_grid[current].y = ( phi[prev_y] - phi[next_y] ) * norm / ( 2.0 * CELL_SIZE );
        }
    }

    // 5. Interpolate back to particles
    mesh_forces.assign( particles.size(), { 0.0, 0.0 } );
    for( size_t i = 0; i < particles.size(); ++i ) {
        const auto& p = particles[i];
        const auto& cd = cic_data[i];

        Vec2 f;
        f.x = force_grid[( cd.ix % MESH_SIZE ) * MESH_SIZE + ( cd.iy % MESH_SIZE )].x * cd.w1 +
            force_grid[( ( cd.ix + 1 ) % MESH_SIZE ) * MESH_SIZE + ( cd.iy % MESH_SIZE )].x * cd.w2 +
            force_grid[( cd.ix % MESH_SIZE ) * MESH_SIZE + ( ( cd.iy + 1 ) % MESH_SIZE )].x * cd.w3 +
            force_grid[( ( cd.ix + 1 ) % MESH_SIZE ) * MESH_SIZE + ( ( cd.iy + 1 ) % MESH_SIZE )].x * cd.w4;

        f.y = force_grid[( cd.ix % MESH_SIZE ) * MESH_SIZE + ( cd.iy % MESH_SIZE )].y * cd.w1 +
            force_grid[( ( cd.ix + 1 ) % MESH_SIZE ) * MESH_SIZE + ( cd.iy % MESH_SIZE )].y * cd.w2 +
            force_grid[( cd.ix % MESH_SIZE ) * MESH_SIZE + ( ( cd.iy + 1 ) % MESH_SIZE )].y * cd.w3 +
            force_grid[( ( cd.ix + 1 ) % MESH_SIZE ) * MESH_SIZE + ( ( cd.iy + 1 ) % MESH_SIZE )].y * cd.w4;

        mesh_forces[i].x = f.x * p.mass;
        mesh_forces[i].y = f.y * p.mass;
    }
}

void compute_PP_forces( std::vector<Particle>& particles, std::vector<Vec2>& pp_forces ) {
    pp_forces.assign( particles.size(), { 0.0, 0.0 } );
    for( size_t i = 0; i < particles.size(); ++i ) {
        for( size_t j = i + 1; j < particles.size(); ++j ) {
            auto& p1 = particles[i];
            auto& p2 = particles[j];

            double dx = displacement( p2.pos.x - p1.pos.x );
            double dy = displacement( p2.pos.y - p1.pos.y );
            double dist_sq = dx * dx + dy * dy;

            if( USE_PM && dist_sq > CUTOFF_RADIUS_SQUARED ) continue;

            double S = 1.0;
            if( USE_PM && dist_sq > R_SWITCH_START_SQ ) {
                double dist = sqrt( dist_sq );
                double x = ( dist - R_SWITCH_START ) / CUTOFF_TRANSITION_WIDTH;
                S = 2 * pow( x, 3 ) - 3 * pow( x, 2 ) + 1;
            }

            double soft_dist_sq = dist_sq + pow( 0.5 * CELL_SIZE, 2 );
            double f_pm_short = G * p1.mass * p2.mass / soft_dist_sq;
            double soft_dist = sqrt( soft_dist_sq );
            Vec2 f_pm_short_vec = { f_pm_short * dx / soft_dist, f_pm_short * dy / soft_dist };

            double pp_dist_sq = dist_sq + SOFTENING_SQUARED;
            double f_pp = G * p1.mass * p2.mass / pp_dist_sq;
            double pp_dist = sqrt( pp_dist_sq );
            Vec2 f_pp_vec = { f_pp * dx / pp_dist, f_pp * dy / pp_dist };

            Vec2 correction_f;
            if( !USE_PM ) {
                f_pm_short_vec.x = 0;
                f_pm_short_vec.y = 0;
            }
            correction_f.x = S * ( f_pp_vec.x - f_pm_short_vec.x );
            correction_f.y = S * ( f_pp_vec.y - f_pm_short_vec.y );

            pp_forces[i].x += correction_f.x;
            pp_forces[i].y += correction_f.y;
            pp_forces[j].x -= correction_f.x;
            pp_forces[j].y -= correction_f.y;
        }
    }
}


void create_zeldovich_ics( std::vector<Particle>& particles ) {
    int n_per_side = N_PER_SIDE;
    int ic_mesh_size = n_per_side;
    double cell_size = DOMAIN_SIZE / ic_mesh_size;

    std::vector<double> real_space_random_field( ic_mesh_size * ic_mesh_size );
    std::default_random_engine generator( SEED );
    std::normal_distribution<double> distribution( 0.0, 1.0 );
    for( auto& val : real_space_random_field ) {
        val = distribution( generator );
    }

    pocketfft::shape_t shape_ic = { ( size_t )ic_mesh_size, ( size_t )ic_mesh_size };
    pocketfft::stride_t stride_r_ic = { static_cast <ptrdiff_t>( sizeof( double ) * ic_mesh_size ), sizeof( double ) };
    pocketfft::stride_t stride_c_ic = { static_cast <ptrdiff_t>( sizeof( std::complex<double> ) ) * ( ic_mesh_size / 2 + 1 ), sizeof( std::complex<double> ) };

    std::vector<std::complex<double>> random_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );
    pocketfft::r2c( shape_ic, stride_r_ic, stride_c_ic, { 0,1 }, true, real_space_random_field.data(), random_k.data(), 1.0 );

    std::vector<std::complex<double>> disp_x_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );
    std::vector<std::complex<double>> disp_y_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );

    double power_spectrum_index = -2.0;

    for( int i = 0; i < ic_mesh_size; ++i ) {
        for( int j = 0; j < ic_mesh_size / 2 + 1; ++j ) {
            double kx_freq = ( i < ic_mesh_size / 2 ) ? i : ( i - ic_mesh_size );
            double ky_freq = j;

            double kx = kx_freq * 2.0 * M_PI / DOMAIN_SIZE;
            double ky = ky_freq * 2.0 * M_PI / DOMAIN_SIZE;
            double k2 = kx * kx + ky * ky;

            int idx = i * ( ic_mesh_size / 2 + 1 ) + j;
            if( i == 0 && j == 0 ) {
                disp_x_k[idx] = { 0,0 };
                disp_y_k[idx] = { 0,0 };
                continue;
            }

            double power_spectrum_sqrt = sqrt( pow( k2, power_spectrum_index / 2.0 ) );
            std::complex<double> delta_k = random_k[idx] * power_spectrum_sqrt;
            std::complex<double> phi_k = -delta_k / k2;

            disp_x_k[idx] = std::complex<double>( 0, 1 ) * kx * phi_k;
            disp_y_k[idx] = std::complex<double>( 0, 1 ) * ky * phi_k;
        }
    }

    std::vector<double> disp_x_real( ic_mesh_size * ic_mesh_size );
    std::vector<double> disp_y_real( ic_mesh_size * ic_mesh_size );

    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, {0,1}, false, disp_x_k.data(), disp_x_real.data(), 1.0);
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, {0,1}, false, disp_y_k.data(), disp_y_real.data(), 1.0 );

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

    double initial_scale_factor = EXPANDING_UNIVERSE ? pow( EXPANSION_START_T, 2.0 / 3.0 ) : 0.5;
    double initial_hubble_param = EXPANDING_UNIVERSE ? ( 2.0 / 3.0 ) / EXPANSION_START_T : 0.0;

    particles.clear();
    double spacing = DOMAIN_SIZE / n_per_side;
    for( int i = 0; i < n_per_side; ++i ) {
        for( int j = 0; j < n_per_side; ++j ) {
            double qx = ( i + 0.5 ) * spacing;
            double qy = ( j + 0.5 ) * spacing;

            double dx = ( disp_x_real[i * n_per_side + j] / std_x ) * initial_scale_factor * cell_size;
            double dy = ( disp_y_real[i * n_per_side + j] / std_y ) * initial_scale_factor * cell_size;

            Particle p;
            p.pos.x = fmod( qx + dx + DOMAIN_SIZE, DOMAIN_SIZE );
            p.pos.y = fmod( qy + dy + DOMAIN_SIZE, DOMAIN_SIZE );
            p.vel.x = initial_hubble_param * dx;
            p.vel.y = initial_hubble_param * dy;
            p.mass = PARTICLE_MASS;
            particles.push_back( p );
        }
    }
}


void KDK_step( double total_time, std::vector<Particle>& particles ) {
    update_cosmology( total_time );
    double a = scale_factor;
    double H = hubble_param;

    for( auto& p : particles ) {
        double total_ax = ( p.acc.x / ( a * a * a ) ) - ( 2 * H * p.vel.x );
        double total_ay = ( p.acc.y / ( a * a * a ) ) - ( 2 * H * p.vel.y );
        p.vel.x += total_ax * DT / 2.0;
        p.vel.y += total_ay * DT / 2.0;
    }

    for( auto& p : particles ) {
        p.pos.x = fmod( p.pos.x + p.vel.x * DT + DOMAIN_SIZE, DOMAIN_SIZE );
        p.pos.y = fmod( p.pos.y + p.vel.y * DT + DOMAIN_SIZE, DOMAIN_SIZE );
    }

    update_cosmology( total_time + DT );
    a = scale_factor;
    H = hubble_param;

    std::vector<Vec2> pp_forces, pm_forces;
    compute_PP_forces( particles, pp_forces );
    compute_mesh_forces( particles, pm_forces );

    for( size_t i = 0; i < particles.size(); ++i ) {
        if( !USE_PM ) {
            pm_forces[i].x = 0;
            pm_forces[i].y = 0;
        }

        auto& p = particles[i];
        Vec2 f = { pp_forces[i].x + pm_forces[i].x, pp_forces[i].y + pm_forces[i].y };
        p.acc.x = f.x / p.mass;
        p.acc.y = f.y / p.mass;

        double total_ax = ( p.acc.x / ( a * a * a ) ) - ( 2 * H * p.vel.x );
        double total_ay = ( p.acc.y / ( a * a * a ) ) - ( 2 * H * p.vel.y );
        p.vel.x += total_ax * DT / 2.0;
        p.vel.y += total_ay * DT / 2.0;
    }
}

double calculate_total_energy( const std::vector<Particle>& particles, double a ) {
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;

    // KE uses proper velocity: v_prop = a * v_peculiar
    for( const auto& p : particles ) {
        double proper_vel_sq = ( a * p.vel.x ) * ( a * p.vel.x ) + ( a * p.vel.y ) * ( a * p.vel.y );
        kinetic_energy += 0.5 * p.mass * proper_vel_sq;
    }

    // PE uses proper distance: r_prop = a * r_comoving
    for( size_t i = 0; i < particles.size(); ++i ) {
        for( size_t j = i + 1; j < particles.size(); ++j ) {
            const auto& p1 = particles[i];
            const auto& p2 = particles[j];

            double dx = displacement( p2.pos.x - p1.pos.x );
            double dy = displacement( p2.pos.y - p1.pos.y );

            // Note: We use the PP softened distance for consistency in measurement
            double proper_dist_sq = ( a * a ) * ( dx * dx + dy * dy + SOFTENING_SQUARED );
            if( proper_dist_sq > 0 ) {
                potential_energy -= G * p1.mass * p2.mass / sqrt( proper_dist_sq );
            }
        }
    }
    return kinetic_energy + potential_energy;
}

Vec2 calculate_total_momentum( const std::vector<Particle>& particles ) {
    Vec2 total_momentum = { 0.0, 0.0 };
    for( const auto& p : particles ) {
        total_momentum.x += p.mass * p.vel.x;
        total_momentum.y += p.mass * p.vel.y;
    }
    return total_momentum;
}

void save_frame( sf::RenderWindow &window, int frame_num ) {
    char filename[100];
    snprintf( filename, sizeof( filename ), "frame_%03d.png", frame_num );

    sf::Texture texture( window.getSize() );
    texture.update( window );
    if( texture.copyToImage().saveToFile( filename ) ) {
        std::cout << "Saved " << filename << std::endl;
    }
}


int main() {
    // --- SFML Window Setup ---
    sf::RenderWindow window( sf::VideoMode( { RENDER_SIZE, RENDER_SIZE } ), "N-Body Simulation" );

    std::vector<Particle> particles;
    create_zeldovich_ics( particles );

    std::vector<Vec2> pp_forces, pm_forces;
    compute_PP_forces( particles, pp_forces );
    compute_mesh_forces( particles, pm_forces );
    for( size_t i = 0; i < particles.size(); ++i ) {
        particles[i].acc.x = ( pp_forces[i].x + pm_forces[i].x ) / particles[i].mass;
        particles[i].acc.y = ( pp_forces[i].y + pm_forces[i].y ) / particles[i].mass;
    }

    int cycle_count = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    update_cosmology( 0 );
    double initial_energy = calculate_total_energy( particles, scale_factor );

    std::cout << "Starting simulation..." << std::endl;
    while( window.isOpen() ) {
        // Handle window events
        while( auto event = window.pollEvent() ) {
            if( event->is<sf::Event::Closed>() ) {
                window.close();
            }
        }

        double total_simulation_time = DT * cycle_count;

        KDK_step( total_simulation_time, particles );

        // --- Rendering ---
        window.clear( sf::Color::Black );

        sf::CircleShape particle_shape( 1.0f );
        particle_shape.setFillColor( sf::Color::White );

        for( const auto& p : particles ) {
            particle_shape.setPosition( { static_cast<float>(p.pos.x * RENDER_SCALE), static_cast<float>(p.pos.y * RENDER_SCALE) } );
            window.draw( particle_shape );
        }

        window.display();

        if( cycle_count % DEBUG_INFO_EVERY_CYCLES == 0 ) {
            update_cosmology( total_simulation_time );
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed_seconds = std::chrono::duration_cast< std::chrono::duration<double> >( now - start_time ).count();
            double energy = calculate_total_energy( particles, scale_factor );
            double energy_error = 100.0 * fabs( energy - initial_energy ) / fabs( initial_energy );

            std::cout << "Cycles: " << cycle_count
                << " | SimTime: " << static_cast< int >( 1000.0 * total_simulation_time )
                << " | ScaleFactor: " << scale_factor
                << " | Energy: " << static_cast<int>( energy_error ) << "%"
                << " | Elapsed Wall Time: " << static_cast< int >( elapsed_seconds ) << "s"
                << std::endl;
        }

        // Save a frame every N cycles
        if( SAVE_RENDER_EVERY_CYCLES != 0 && cycle_count > 0 && cycle_count % SAVE_RENDER_EVERY_CYCLES == 0 ) {
            int frame_num = cycle_count / SAVE_RENDER_EVERY_CYCLES;
            save_frame( window, frame_num );
        }

        cycle_count++;
    }

    return 0;
}

