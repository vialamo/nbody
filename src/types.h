#pragma once
#include <Eigen/Dense>

namespace TimeConstants {
// The smallest valid macro-step the engine can take.
// Steps smaller than this are considered rounding noise
constexpr double MIN_VALID_DT = 1e-10;

// The floating-point precision limit for checking if two times are aligned.
// Used to check if a particle's local clock matches the global clock
constexpr double SYNC_TOLERANCE = 1e-12;

// The tolerance for verifying if a continuous global time matches
// an integer multiple of a fixed time bin or snapshot interval
constexpr double BIN_BOUNDARY_TOLERANCE = 1e-7;

constexpr int MAX_TIMESTEP_BIN_DEPTH = 30;
}  // namespace TimeConstants

struct Vec3 {
    double x, y, z;
};

struct Vec3i {
    int x, y, z;
};

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

struct CIC_Data {
    int ix, iy, iz;
    // 8 weights for the 8 corners of a 3D cube
    double w000, w100, w010, w110, w001, w101, w011, w111;
};

struct TimestepInfo {
    double dt_macro;
    double dt_hydro;
    double dt_grav;
    double dt_cool;
};

struct EnergyBudget {
    double radiated = 0.0;
    double photoheated = 0.0;
    double gravitational_work = 0.0;
    double expansion_work = 0.0;
    double entropy_switch = 0.0;
    double floor_injection = 0.0;

    // Allow addition of deltas
    EnergyBudget& operator+=(const EnergyBudget& rhs) {
        radiated += rhs.radiated;
        photoheated += rhs.photoheated;
        gravitational_work += rhs.gravitational_work;
        expansion_work += rhs.expansion_work;
        entropy_switch += rhs.entropy_switch;
        floor_injection += rhs.floor_injection;
        return *this;
    }
};
