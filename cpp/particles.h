#pragma once
#include "types.h"
#include "config.h"
#include <vector>
#include <utility>

class ParticleSystem {
public:
    std::vector<Particle> particles;

    Grid dm_rho;
    std::vector<CIC_Data> cic_data;
    CellGrid cell_grid;

    ParticleSystem( const Config& config );

    void bin_and_assign_mass( const Config& config );

    void interpolate_cic_forces( const Grid& ax_grid, const Grid& ay_grid,
        std::vector<Vec2>& forces, const Config& config );

    void compute_pp_forces( std::vector<Vec2>& pp_forces, const Config& config );

    std::pair<double, double> calculate_energies( double a, const Config& config ) const;

    double get_gravity_timestep( const Config& config ) const;
};