#pragma once
#include "types.h"
#include "config.h"
#include <vector>
#include <utility>

class ParticleSystem {
public:
    std::vector<Particle> particles;

    Grid3D dm_rho;
    std::vector<CIC_Data> cic_data;
    CellGrid cell_grid;

    ParticleSystem( const Config& config );

    void bin_and_assign_mass( const Config& config );

    void interpolate_cic_forces( const Grid3D& ax_grid, const Grid3D& ay_grid, const Grid3D& az_grid,
        std::vector<Vec3>& forces, const Config& config );

    void compute_pp_forces( std::vector<Vec3>& pp_forces, const Config& config );

    double calculate_kinetic_energy(double a) const;

    double get_gravity_timestep( const Config& config ) const;
};