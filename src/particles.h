#pragma once
#include <vector>

#include "config.h"
#include "types.h"

// Encapsulates a flat-array linked list for particle binning
struct CellLinkedList {
    // Size = num_cells. Index of the first particle in a given cell
    std::vector<int> first_in_cell;

    // Size = num_particles. Index of the next particle sharing the same cell
    std::vector<int> next_particle;

    void resize(int num_cells, int num_particles) {
        first_in_cell.assign(num_cells, -1);
        next_particle.assign(num_particles, -1);
    }

    void reset() {
        std::fill(first_in_cell.begin(), first_in_cell.end(), -1);
        std::fill(next_particle.begin(), next_particle.end(), -1);
    }
};

class ParticleSystem {
   private:
    CellLinkedList cell_list;

    std::vector<CIC_Data> cic_data;

    friend struct ParticleTestAccess;

   public:
    Grid3D dm_rho;

    size_t num_particles = 0;

    std::vector<double> pos_x;
    std::vector<double> pos_y;
    std::vector<double> pos_z;

    std::vector<double> vel_x;
    std::vector<double> vel_y;
    std::vector<double> vel_z;

    std::vector<double> acc_x;
    std::vector<double> acc_y;
    std::vector<double> acc_z;

    std::vector<double> mass;

    double max_accel_sq;

    ParticleSystem(const Config& config);

    void bin_and_assign_mass(const Config& config);

    void interpolate_cic_forces(const Grid3D& ax_grid, const Grid3D& ay_grid,
                                const Grid3D& az_grid, const Config& config);

    void compute_pp_forces(const Config& config);

    double calculate_kinetic_energy(double a) const;

    double get_gravity_timestep(const Config& config) const;

    const Grid3D& get_rho() const { return dm_rho; }

    void add_particle(double px, double py, double pz, double vx, double vy,
                      double vz, double m);
};