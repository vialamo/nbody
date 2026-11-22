#pragma once

#include <cstdint>

namespace TimeConstants {
constexpr int MAX_TIMESTEP_BIN_DEPTH = 30;
constexpr uint64_t TIMEBASE = 1ULL << MAX_TIMESTEP_BIN_DEPTH;
constexpr int MAX_NGB_TIME_FACTOR = 4;
}  // namespace TimeConstants

struct Vec3 {
    double x, y, z;
};

struct Vec3i {
    int x, y, z;
};

struct CIC_Data {
    int ix, iy, iz;
    // 8 weights for the 8 corners of a 3D cube
    double w000, w100, w010, w110, w001, w101, w011, w111;
};

struct TimestepInfo {
    uint64_t ti_macro;
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
    double def_switch = 0.0;
    double entropy_switch = 0.0;
    double floor_injection = 0.0;

    // Allow addition of deltas
    EnergyBudget& operator+=(const EnergyBudget& rhs) {
        radiated += rhs.radiated;
        photoheated += rhs.photoheated;
        gravitational_work += rhs.gravitational_work;
        expansion_work += rhs.expansion_work;
        def_switch += rhs.def_switch;
        entropy_switch += rhs.entropy_switch;
        floor_injection += rhs.floor_injection;
        return *this;
    }
};
