#pragma once

#include <fftw3.h>

#include <Eigen/Dense>
#include <vector>

#include "types.h"

struct Config;
struct SimState;

class PM {
   private:
    std::vector<std::complex<double>> rho_k;
    std::vector<std::complex<double>> phi_k;
    fftw_plan forward_plan = nullptr;
    fftw_plan backward_plan = nullptr;

   public:
    struct Grid3D {
        int nx, ny, nz;
        Eigen::VectorXd data;

        Grid3D() : nx(0), ny(0), nz(0) {}

        Grid3D(const Vec3i& size)
            : nx(size.x), ny(size.y), nz(size.z), data(nx * ny * nz) {
            data.setZero();
        }

        // Access operator (row-major mapping: x varies slowest, z fastest)
        inline double& operator()(int x, int y, int z) {
            return data[x * ny * nz + y * nz + z];
        }

        inline const double& operator()(int x, int y, int z) const {
            return data[x * ny * nz + y * nz + z];
        }

        // Pass-throughs for Eigen array operations (Unchanged)
        auto array() { return data.array(); }
        auto array() const { return data.array(); }
        void setZero() { data.setZero(); }
        void fill(double val) { data.setConstant(val); }
        double minCoeff() const { return data.minCoeff(); }
        double maxCoeff() const { return data.maxCoeff(); }
        double sum() const { return data.sum(); }

        Grid3D roll(int shift, int axis) const;

        inline size_t size() const { return data.size(); }

        // For HDF5 pointers
        const double* raw_data() const { return data.data(); }
    };

    Grid3D rho;
    Grid3D phi;

    // Comoving gravitational acceleration on the mesh
    Grid3D acc_x;
    Grid3D acc_y;
    Grid3D acc_z;

    PM(const Config& config);
    ~PM();

    // Splats particle masses onto a grid and caches the interpolation weights
    void bin_and_assign_mass(SimState& state, const Config& config);

    // Takes the density grid, solves Poisson's equation in Fourier space,
    // and computes the spatial gradients to yield acceleration grids
    void compute_acceleration(SimState& state, const Config& config);

    // Gathers forces from grids back to the particles using cached weights
    void interpolate_forces(SimState& state, const Config& config);
};
