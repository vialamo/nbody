#pragma once
#include <Eigen/Dense>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Calculates the shortest distance between two points in a periodic domain
#ifdef USE_GPU
#pragma omp declare target
#endif
inline double periodic_displacement(double dx, double domain_size) {
    double half_domain = 0.5 * domain_size;
    dx = (dx > half_domain) ? dx - domain_size : dx;
    dx = (dx < -half_domain) ? dx + domain_size : dx;
    return dx;
}
#ifdef USE_GPU
#pragma omp end declare target
#endif

class QuadraticInterpolator {
   private:
    double a, b, c;

   public:
    QuadraticInterpolator(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                          const Eigen::Vector2d& p3);

    inline double evaluate(double x) const { return a * x * x + b * x + c; }
};

// Expands a 21-bit integer into 64 bits by inserting 2 zeros after each bit.
inline uint64_t splitBy3(uint32_t a) {
    uint64_t x = a & 0x1fffff;  // 21 bits
    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8) & 0x100f00f00f00f00f;
    x = (x | x << 4) & 0x10c30c30c30c30c3;
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}

// Combines three 21-bit coordinates into a single 63-bit Morton code
inline uint64_t morton3D(uint32_t x, uint32_t y, uint32_t z) {
    return (splitBy3(x) | (splitBy3(y) << 1) | (splitBy3(z) << 2));
}

#pragma omp declare target
inline double min_periodic_dist_sq(double p, double min_b, double max_b,
                                   double domain) {
    // Standard distance to AABB
    double d = 0.0;
    if (p < min_b)
        d = min_b - p;
    else if (p > max_b)
        d = p - max_b;

    // Check wrapped distances
    double d_wrapped_left = 0.0;
    double p_wrapped_left = p + domain;
    if (p_wrapped_left < min_b)
        d_wrapped_left = min_b - p_wrapped_left;
    else if (p_wrapped_left > max_b)
        d_wrapped_left = p_wrapped_left - max_b;

    double d_wrapped_right = 0.0;
    double p_wrapped_right = p - domain;
    if (p_wrapped_right < min_b)
        d_wrapped_right = min_b - p_wrapped_right;
    else if (p_wrapped_right > max_b)
        d_wrapped_right = p_wrapped_right - max_b;

    double min_d = std::min({d, d_wrapped_left, d_wrapped_right});
    return min_d * min_d;
}
#pragma omp end declare target