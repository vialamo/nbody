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
    dx = (dx >  half_domain) ? dx - domain_size : dx;
    dx = (dx < -half_domain) ? dx + domain_size : dx;
    return dx;
}
#ifdef USE_GPU
#pragma omp end declare target
#endif