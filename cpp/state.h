#pragma once
#include "gas.h"
#include "particles.h"
#include "types.h"

struct SimState {
    ParticleSystem dm;
    GasGrid gas;

    // Global Gravitational Field on the Eulerian Mesh
    Grid3D gravity_x;
    Grid3D gravity_y;
    Grid3D gravity_z;

    double total_time = 0.0;
    double scale_factor = 1.0;
    double hubble_param = 0.0;

    SimState(const Config& config)
        : dm(config),
          gas(config),
          gravity_x(config.MESH_SIZE),
          gravity_y(config.MESH_SIZE),
          gravity_z(config.MESH_SIZE) {}
};