import os
import math
import numpy as np
import configparser

class Config:
    """Loads all simulation parameters from a .ini file."""
    def __init__(self, config_filename):
        config = configparser.ConfigParser()
        if not os.path.exists(config_filename):
            raise FileNotFoundError(f"Config file '{config_filename}' not found.")
        
        config.read(config_filename)
        self.params = {} 
        
        # [domain]
        self.domain_size = self._get(config, 'domain', 'domain_size', 'float')
        self.mesh_size = self._get(config, 'domain', 'mesh_size', 'int')
        
        # [matter]
        self.omega_baryon = self._get(config, 'matter', 'omega_baryon', 'float')
        self.num_dm_particles_side = self._get(config, 'matter', 'num_dm_particles_side', 'int')
        self.enable_hydro = self._get(config, 'matter', 'enable_hydro', 'bool')
        
        # [physics]
        self.expansion_start_t = self._get(config, 'physics', 'expansion_start_t', 'float')
        self.expanding_universe = self._get(config, 'physics', 'expanding_universe', 'bool')
        self.g_const = self._get(config, 'physics', 'g_const', 'float') 
        self.gamma = self._get(config, 'physics', 'gamma', 'float')
        self.initial_power_spectrum_index = self._get(config, 'physics', 'initial_power_spectrum_index', 'float')

        # [p3m]
        self.cutoff_radius_cells = self._get(config, 'p3m', 'cutoff_radius_cells', 'float')
        self.cutoff_transition_width_factor = self._get(config, 'p3m', 'cutoff_transition_width_factor', 'float')

        # [time]
        self.dt_factor = self._get(config, 'time', 'dt_factor', 'float')
        self.cfl_safety_factor = self._get(config, 'time', 'cfl_safety_factor', 'float')
        self.gravity_dt_factor = self._get(config, 'time', 'gravity_dt_factor', 'float')
        self.use_adaptive_dt = self._get(config, 'time', 'use_adaptive_dt', 'bool')
        self.max_cycles = self._get(config, 'time', 'max_cycles', 'int', default=1000)

        # [output]
        self.debug_info_every_cycles = self._get(config, 'output', 'debug_info_every_cycles', 'int')
        self.save_snapshot_every_cycles = self._get(config, 'output', 'save_snapshot_every_cycles', 'int')
        self.seed = self._get(config, 'output', 'seed', 'int')

        # [visualization]
        self.render_size = self._get(config, 'visualization', 'render_size', 'int')
        self.enable_visualization = self._get(config, 'visualization', 'enable_visualization', 'bool', default=True)

        # --- Calculate derived parameters ---
        self.cell_size = self.domain_size / self.mesh_size
        self.cell_volume = self.cell_size**3
        self.omega_dm = 1.0 - self.omega_baryon
        
        self.cutoff_radius = self.cutoff_radius_cells * self.cell_size
        self.cutoff_radius_squared = self.cutoff_radius**2
        self.cutoff_transition_width = self.cutoff_transition_width_factor * self.cutoff_radius
        self.r_switch_start = self.cutoff_radius - self.cutoff_transition_width
        self.r_switch_start_sq = self.r_switch_start**2
        
        self.num_dm_particles = self.num_dm_particles_side**3
        self.dm_particle_mass = self.omega_dm / self.num_dm_particles
        self.gas_total_mass = self.omega_baryon
        
        self.mean_interparticle_spacing = self.domain_size / (self.num_dm_particles**(1/3))
        self.softening_squared = (self.mean_interparticle_spacing / 50.0)**2
        
        self.initial_hubble_param = (2.0/3.0) / self.expansion_start_t if self.expanding_universe else 0.0
        self.initial_scale_factor = self.expansion_start_t**(2.0/3.0) if self.expanding_universe else 0.0
        
        self.dynamical_time = 1.0 / math.sqrt(4.0 * math.pi * self.g_const)
        self.fixed_dt = self.dt_factor * self.dynamical_time
        
        self.render_scale = self.render_size / self.domain_size

        np.random.seed(self.seed)

    def _get(self, config, section, option, dtype='str', default=None):
        try:
            if dtype == 'float':
                val = config.getfloat(section, option)
            elif dtype == 'int':
                val = config.getint(section, option)
            elif dtype == 'bool':
                val = config.getboolean(section, option)
            else:
                val = config.get(section, option)
        except (configparser.NoOptionError, configparser.NoSectionError):
            if default is not None:
                val = default
            else:
                raise
        self.params[f"{section}_{option}"] = val
        return val

    def get_all_params(self):
        return self.params