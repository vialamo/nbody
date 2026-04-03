#pragma once
#include <vector>
#include <Eigen/Dense>

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

struct Particle {
    Vec3 pos;
    Vec3 vel;
    Vec3 acc;
    double mass;
};

struct Grid3D {
    int n;
    Eigen::VectorXd data;

    // Default constructor needed for structs/arrays
    Grid3D() : n( 0 ) {}

    Grid3D( int size ) : n( size ), data( size* size* size ) {
        data.setZero();
    }

    // Access operator (row-major mapping)
    inline double& operator()( int x, int y, int z ) {
        return data[x * n * n + y * n + z];
    }

    inline const double& operator()( int x, int y, int z ) const {
        return data[x * n * n + y * n + z];
    }

    // Pass-throughs for Eigen array operations
    auto array() { return data.array(); }
    auto array() const { return data.array(); }
    void setZero() { data.setZero(); }
    void fill( double val ) { data.setConstant( val ); }
    double minCoeff() const { return data.minCoeff(); }
    double maxCoeff() const { return data.maxCoeff(); }
    double sum() const { return data.sum(); }

    // For HDF5 pointers
    const double* raw_data() const { return data.data(); }
};

// Use a flat std::vector of vectors for the 3D cell grid
using CellGrid = std::vector<std::vector<int>>;

struct CIC_Data {
    int ix, iy, iz;
    // 8 weights for the 8 corners of a 3D cube
    double w000, w100, w010, w110, w001, w101, w011, w111;
};