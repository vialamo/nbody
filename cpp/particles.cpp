#include "particles.h"
#include <cmath>
#include <limits>

static double displacement( double dx, const Config& config ) {
    return fmod( dx + 0.5 * config.DOMAIN_SIZE + config.DOMAIN_SIZE, config.DOMAIN_SIZE ) - 0.5 * config.DOMAIN_SIZE;
}

ParticleSystem::ParticleSystem( const Config& config )
    : dm_rho( config.MESH_SIZE, config.MESH_SIZE ),
    cic_data( config.NUM_DM_PARTICLES ),
    cell_grid( config.MESH_SIZE* config.MESH_SIZE )
{
    // Pre-allocate memory
    particles.reserve( config.NUM_DM_PARTICLES );
}

void ParticleSystem::bin_and_assign_mass( const Config& config ) {
    dm_rho.setZero();
    cic_data.assign( particles.size(), {} );

    for( auto& cell : cell_grid ) {
        cell.clear();
    }

    for( size_t i = 0; i < particles.size(); ++i ) {
        const auto& p = particles[i];
        int ix = static_cast< int >( p.pos.x / config.CELL_SIZE );
        int iy = static_cast< int >( p.pos.y / config.CELL_SIZE );
        double frac_x = ( p.pos.x / config.CELL_SIZE ) - ix;
        double frac_y = ( p.pos.y / config.CELL_SIZE ) - iy;

        double w1 = ( 1 - frac_x ) * ( 1 - frac_y ), w2 = frac_x * ( 1 - frac_y );
        double w3 = ( 1 - frac_x ) * frac_y, w4 = frac_x * frac_y;
        cic_data[i] = { ix, iy, w1, w2, w3, w4 };

        dm_rho( ( ix + config.MESH_SIZE ) % config.MESH_SIZE, ( iy + config.MESH_SIZE ) % config.MESH_SIZE ) += p.mass * w1;
        dm_rho( ( ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( iy + config.MESH_SIZE ) % config.MESH_SIZE ) += p.mass * w2;
        dm_rho( ( ix + config.MESH_SIZE ) % config.MESH_SIZE, ( iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) += p.mass * w3;
        dm_rho( ( ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) += p.mass * w4;

        int cell_index_x = ( ix % config.MESH_SIZE + config.MESH_SIZE ) % config.MESH_SIZE;
        int cell_index_y = ( iy % config.MESH_SIZE + config.MESH_SIZE ) % config.MESH_SIZE;
        int cell_index = cell_index_y * config.MESH_SIZE + cell_index_x;
        cell_grid[cell_index].push_back( static_cast< int >( i ) );
    }

    dm_rho /= ( config.CELL_SIZE * config.CELL_SIZE );
}

void ParticleSystem::interpolate_cic_forces( const Grid& ax_grid, const Grid& ay_grid, std::vector<Vec2>& forces, const Config& config ) {
    forces.assign( particles.size(), { 0.0, 0.0 } );
    for( size_t i = 0; i < particles.size(); ++i ) {
        const auto& p = particles[i];
        const auto& cd = cic_data[i];
        double ax = ( ax_grid( ( cd.ix + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w1 +
            ax_grid( ( cd.ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w2 +
            ax_grid( ( cd.ix + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w3 +
            ax_grid( ( cd.ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w4 );
        double ay = ( ay_grid( ( cd.ix + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w1 +
            ay_grid( ( cd.ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w2 +
            ay_grid( ( cd.ix + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w3 +
            ay_grid( ( cd.ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w4 );
        forces[i].x = ax * p.mass;
        forces[i].y = ay * p.mass;
    }
}

void ParticleSystem::compute_pp_forces( std::vector<Vec2>& pp_forces, const Config& config ) {
    pp_forces.assign( particles.size(), { 0.0, 0.0 } );
    const int search_radius_cells = static_cast< int >( ceil( config.CUTOFF_RADIUS / config.CELL_SIZE ) );

    for( size_t i = 0; i < particles.size(); ++i ) {
        auto& p1 = particles[i];
        int ix = static_cast< int >( p1.pos.x / config.CELL_SIZE );
        int iy = static_cast< int >( p1.pos.y / config.CELL_SIZE );

        for( int dx_cell = -search_radius_cells; dx_cell <= search_radius_cells; ++dx_cell ) {
            for( int dy_cell = -search_radius_cells; dy_cell <= search_radius_cells; ++dy_cell ) {
                int neighbor_ix = ( ix + dx_cell + config.MESH_SIZE ) % config.MESH_SIZE;
                int neighbor_iy = ( iy + dy_cell + config.MESH_SIZE ) % config.MESH_SIZE;
                int neighbor_cell_index = neighbor_iy * config.MESH_SIZE + neighbor_ix;

                for( int j : cell_grid[neighbor_cell_index] ) {
                    if( i >= j ) continue;
                    auto& p2 = particles[j];

                    double dx = displacement( p2.pos.x - p1.pos.x, config );
                    double dy = displacement( p2.pos.y - p1.pos.y, config );
                    double dist_sq = dx * dx + dy * dy;

                    if( config.USE_PM && dist_sq > config.CUTOFF_RADIUS_SQUARED ) continue;

                    double S = 1.0;
                    if( config.USE_PM && dist_sq > config.R_SWITCH_START_SQ ) {
                        double dist = sqrt( dist_sq );
                        double x = ( dist - config.R_SWITCH_START ) / config.CUTOFF_TRANSITION_WIDTH;
                        S = 2 * pow( x, 3 ) - 3 * pow( x, 2 ) + 1;
                    }

                    double soft_dist_sq = dist_sq + pow( 0.5 * config.CELL_SIZE, 2 );
                    double f_pm_short = config.G * p1.mass * p2.mass / soft_dist_sq;
                    double soft_dist = sqrt( soft_dist_sq );
                    Vec2 f_pm_short_vec = { f_pm_short * dx / soft_dist, f_pm_short * dy / soft_dist };

                    double pp_dist_sq = dist_sq + config.SOFTENING_SQUARED;
                    double f_pp = config.G * p1.mass * p2.mass / pp_dist_sq;
                    double pp_dist = sqrt( pp_dist_sq );
                    Vec2 f_pp_vec = { f_pp * dx / pp_dist, f_pp * dy / pp_dist };

                    Vec2 correction_f;
                    if( !config.USE_PM ) {
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
    }
}

std::pair<double, double> ParticleSystem::calculate_energies( double a, const Config& config ) const {
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;
    for( const auto& p : particles ) {
        double proper_vel_sq = ( a * p.vel.x ) * ( a * p.vel.x ) + ( a * p.vel.y ) * ( a * p.vel.y );
        kinetic_energy += 0.5 * p.mass * proper_vel_sq;
    }
    for( size_t i = 0; i < particles.size(); ++i ) {
        for( size_t j = i + 1; j < particles.size(); ++j ) {
            const auto& p1 = particles[i];
            const auto& p2 = particles[j];
            double dx = displacement( p2.pos.x - p1.pos.x, config );
            double dy = displacement( p2.pos.y - p1.pos.y, config );
            double proper_dist_sq = ( a * a ) * ( dx * dx + dy * dy + config.SOFTENING_SQUARED );
            if( proper_dist_sq > 0 ) {
                potential_energy -= config.G * p1.mass * p2.mass / sqrt( proper_dist_sq );
            }
        }
    }
    return { kinetic_energy, potential_energy };
}

double ParticleSystem::get_gravity_timestep( const Config& config ) const {
    if( particles.empty() ) return std::numeric_limits<double>::infinity();

    double max_accel_sq = 1e-9;
    for( const auto& p : particles ) {
        double accel_sq = p.acc.x * p.acc.x + p.acc.y * p.acc.y;
        if( accel_sq > max_accel_sq ) {
            max_accel_sq = accel_sq;
        }
    }
    double dt_grav = sqrt( config.SOFTENING_SQUARED ) / sqrt( max_accel_sq );
    return dt_grav * config.GRAVITY_DT_FACTOR;
}