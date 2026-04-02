#pragma once
#include <vector>
#include <Eigen/Dense>

// Global Eigen typedef
using Grid = Eigen::MatrixXd;

struct Vec2 {
    double x = 0.0, y = 0.0;
};

struct Particle {
    Vec2 pos;
    Vec2 vel;
    Vec2 acc;
    double mass;
};

// Data structure for the PP solver's cell list
using CellGrid = std::vector<std::vector<int>>;

struct CIC_Data {
    int ix, iy;
    double w1, w2, w3, w4;
};