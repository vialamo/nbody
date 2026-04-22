#include <cmath>

// Calculates the shortest distance between two points in a periodic domain
inline double periodic_displacement(double dx, double domain_size) {
    return std::fmod(dx + 0.5 * domain_size + domain_size, domain_size) - 0.5 * domain_size;
}