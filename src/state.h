#pragma once
#include <memory>

#include "cooling.h"
#include "gas.h"
#include "mfm.h"
#include "particles.h"
#include "types.h"

struct SimState {
    ParticleSystem dm;
    std::unique_ptr<GasGrid> gas;
    std::unique_ptr<GasParticleSystem> mfm_gas;
    Cooling cooling;

    Grid3D total_rho;
    Grid3D phi;

    // Global Gravitational Field on the Eulerian Mesh
    Grid3D pm_gravity_x;
    Grid3D pm_gravity_y;
    Grid3D pm_gravity_z;

    double total_time = 0.0;
    double scale_factor = 1.0;
    double hubble_param = 0.0;

    SimState(const Config& config)
        : dm(config),
          gas(nullptr),
          mfm_gas(nullptr),
          cooling(config),
          pm_gravity_x(config.mesh_size),
          pm_gravity_y(config.mesh_size),
          pm_gravity_z(config.mesh_size),
          total_rho(config.mesh_size),
          phi(config.mesh_size) {
        if (config.hydro_method == HydroMethod::Eulerian) {
            gas = std::make_unique<GasGrid>(config);
        } else if (config.hydro_method == HydroMethod::MFM) {
            mfm_gas = std::make_unique<GasParticleSystem>(config);
        }
    }
};