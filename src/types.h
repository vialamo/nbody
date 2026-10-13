#pragma once
#include <Eigen/Dense>
#include <vector>

struct Vec3 {
    double x, y, z;
};

struct Grid3D {
    int nx, ny, nz;
    Eigen::VectorXd data;

    Grid3D() : nx(0), ny(0), nz(0) {}

    Grid3D(const Vec3& size)
        : nx(static_cast<int>(size.x)),
          ny(static_cast<int>(size.y)),
          nz(static_cast<int>(size.z)),
          data(nx * ny * nz) {
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
