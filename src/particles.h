#pragma once
#include <vector>

#include "config.h"
#include "types.h"

class ParticleSystem {
   private:
    CellGrid cell_grid;
    std::vector<CIC_Data> cic_data;
    Grid3D dm_rho;

    friend struct ParticleTestAccess;

   public:
    std::vector<Particle> particles;
    double max_accel_sq;

    ParticleSystem(const Config& config);

    void bin_and_assign_mass(const Config& config);

    void interpolate_cic_forces(const Grid3D& ax_grid, const Grid3D& ay_grid,
                                const Grid3D& az_grid,
                                std::vector<Vec3>& forces,
                                const Config& config);

    void compute_pp_forces(std::vector<Vec3>& pp_forces, const Config& config);

    double calculate_kinetic_energy(double a) const;

    double get_gravity_timestep(const Config& config) const;

    const std::vector<Particle>& get_particles() const { return particles; }
    const Grid3D& get_rho() const { return dm_rho; }
};