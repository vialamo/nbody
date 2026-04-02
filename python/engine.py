import time
from ics import initialize_state 
from integrator import KDK_step, calculate_diagnostics
from gas import get_cfl_timestep
from particles import get_gravity_timestep
from utils import save_snapshot

class SimulationEngine:
    def __init__(self, config, logger, h5_file):
        self.config = config
        self.logger = logger
        self.h5_file = h5_file
        
        self.cycle_count = 0
        self.snapshot_count = 0
        
        self.state = initialize_state(config)
        
        # Setup initial dt
        self.dt_cfl = get_cfl_timestep(self.state.gas, self.config)
        self.dt_grav = get_gravity_timestep(self.state.particles, self.config)
        if self.config.use_adaptive_dt:
            self.current_dt = min(self.dt_cfl, self.dt_grav, self.config.fixed_dt)
        else:
            self.current_dt = self.config.fixed_dt

        # Log initial state
        diag = calculate_diagnostics(self.state.particles, self.state.gas, {}, self.current_dt, 
                                     0, self.state.total_time, self.state.scale_factor, self.config)
        self.logger.log(diag)

    def step(self):
        timings = KDK_step(self.state, self.current_dt, self.config)
        self.cycle_count += 1
        
        # Update Timestep for next cycle
        if self.config.use_adaptive_dt:
            self.dt_cfl = get_cfl_timestep(self.state.gas, self.config)
            self.dt_grav = get_gravity_timestep(self.state.particles, self.config)
            self.current_dt = min(self.dt_cfl, self.dt_grav, self.config.fixed_dt)
        else:
            self.current_dt = self.config.fixed_dt

        # I/O and Logging (Reading purely from the state object)
        t_io_start = time.time()
        if self.config.save_snapshot_every_cycles > 0 and self.cycle_count % self.config.save_snapshot_every_cycles == 0:
            snapshot_name = f"snapshot_{self.snapshot_count:04d}"
            save_snapshot(self.h5_file, snapshot_name, self.state.particles, self.state.gas, 
                          self.state.total_time, self.state.scale_factor, self.config)
            self.snapshot_count += 1
        timings['io'] = time.time() - t_io_start

        if self.config.debug_info_every_cycles > 0 and self.cycle_count % self.config.debug_info_every_cycles == 0:
            diag = calculate_diagnostics(self.state.particles, self.state.gas, timings, 
                                         self.current_dt, self.cycle_count, 
                                         self.state.total_time, self.state.scale_factor, self.config)
            self.logger.log(diag)