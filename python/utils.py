import time
import csv
import numpy as np

class Diagnostics:
    def __init__(self):
        self.cycle = 0
        self.sim_time = 0.0
        self.scale_factor = 1.0
        self.total_mass_gas = 0.0
        self.total_mass_dm = 0.0
        self.total_momentum_gas = (0.0, 0.0, 0.0)
        self.total_momentum_dm = (0.0, 0.0, 0.0)
        self.ke_gas = 0.0
        self.ke_dm = 0.0
        self.pe_dm = 0.0
        self.ie_gas = 0.0
        self.dt_cfl = 0.0
        self.dt_gravity = 0.0
        self.dt_final = 0.0
        self.max_gas_density = 0.0
        self.max_gas_pressure = 0.0
        self.max_gas_velocity = 0.0
        self.wall_time_total = 0.0
        self.wall_time_pm = 0.0
        self.wall_time_pp = 0.0
        self.wall_time_hydro = 0.0
        self.wall_time_io = 0.0
        
    @property
    def total_mass(self):
        return self.total_mass_gas + self.total_mass_dm

    @property
    def total_momentum(self):
        return (self.total_momentum_gas[0] + self.total_momentum_dm[0],
                self.total_momentum_gas[1] + self.total_momentum_dm[1],
                self.total_momentum_gas[2] + self.total_momentum_dm[2])
    
    @property
    def total_energy(self):
        return self.ke_gas + self.ke_dm + self.pe_dm + self.ie_gas

class Logger:
    def __init__(self, log_filename="diagnostics.csv"):
        self.log_filename = log_filename
        self.log_file = open(log_filename, 'w', newline='')
        self.writer = None
        self.start_time = time.time()

    def write_header(self, diag):
        headers = [
            'cycle', 'sim_time', 'scale_factor',
            'total_mass_gas', 'total_mass_dm',
            'total_momentum_gas_x', 'total_momentum_gas_y', 'total_momentum_gas_z',
            'total_momentum_dm_x', 'total_momentum_dm_y', 'total_momentum_dm_z',
            'ke_gas', 'ke_dm', 'pe_dm', 'ie_gas',
            'dt_cfl', 'dt_gravity', 'dt_final',
            'max_gas_density', 'max_gas_pressure', 'max_gas_velocity',
            'wall_time_total', 'wall_time_pm', 'wall_time_pp', 'wall_time_hydro', 'wall_time_io'
        ]
        flat_diag = self._flatten_diag(diag)
        self.fieldnames = [h for h in headers if h in flat_diag]
        self.writer = csv.DictWriter(self.log_file, fieldnames=self.fieldnames)
        self.writer.writeheader()
        self.log_file.flush()

    def _flatten_diag(self, diag):
        flat = diag.__dict__.copy()
        mom_g = flat.pop('total_momentum_gas')
        flat['total_momentum_gas_x'] = mom_g[0]
        flat['total_momentum_gas_y'] = mom_g[1]
        flat['total_momentum_gas_z'] = mom_g[2]
        mom_dm = flat.pop('total_momentum_dm')
        flat['total_momentum_dm_x'] = mom_dm[0]
        flat['total_momentum_dm_y'] = mom_dm[1]
        flat['total_momentum_dm_z'] = mom_dm[2]
        return flat

    def log(self, diag):
        if self.writer is None:
            self.write_header(diag)
        flat_diag = self._flatten_diag(diag)
        csv_row = {k: v for k, v in flat_diag.items() if k in self.fieldnames}
        self.writer.writerow(csv_row)
        self.log_file.flush()
        
        wall_time_str = f"Wall: {time.time() - self.start_time:.1f}s"
        print(f"[Cycle {diag.cycle}] {wall_time_str} | SimTime: {diag.sim_time:.3f} | a: {diag.scale_factor:.3f}")
        mass_err = (diag.total_mass - 1.0)
        print(f"  [Physics]\n    - Mass (P/G/T):   {diag.total_mass_dm:.4f} | {diag.total_mass_gas:.4f} | {diag.total_mass:.4f} (Err: {mass_err:+.1e})")
        print(f"    - Energy (Total): {diag.total_energy:.3e}")
        print(f"  [Stability]\n    - Timestep (Final): {diag.dt_final:.2e}")
        print("-" * 70)

    def close(self):
        self.log_file.close()

def save_snapshot(h5_file, snapshot_name, particles, gas, sim_time, scale_factor, config):
    try:
        grp = h5_file.create_group(snapshot_name)
        grp.attrs['time'] = sim_time
        grp.attrs['scale_factor'] = scale_factor
        
        num_p = len(particles)
        positions = np.zeros((num_p, 3))
        velocities = np.zeros((num_p, 3))
        masses = np.zeros(num_p)
        for i, p in enumerate(particles):
            positions[i] = [p.x, p.y, p.z]
            velocities[i] = [p.vx, p.vy, p.vz]
            masses[i] = p.mass
            
        p_grp = grp.create_group("particles")
        p_grp.create_dataset("positions", data=positions, compression="gzip")
        p_grp.create_dataset("velocities", data=velocities, compression="gzip")
        p_grp.create_dataset("masses", data=masses, compression="gzip")
        
        if config.enable_hydro:
            g_grp = grp.create_group("gas")
            g_grp.create_dataset("density", data=gas.density, compression="gzip")
            g_grp.create_dataset("momentum_x", data=gas.momentum_x, compression="gzip")
            g_grp.create_dataset("momentum_y", data=gas.momentum_y, compression="gzip")
            g_grp.create_dataset("momentum_z", data=gas.momentum_z, compression="gzip")
            g_grp.create_dataset("energy", data=gas.energy, compression="gzip")
            g_grp.create_dataset("pressure", data=gas.pressure, compression="gzip") 
            g_grp.create_dataset("velocity_x", data=gas.velocity_x, compression="gzip")
            g_grp.create_dataset("velocity_y", data=gas.velocity_y, compression="gzip")
            g_grp.create_dataset("velocity_z", data=gas.velocity_z, compression="gzip")
    except Exception as e:
        print(f"Error saving snapshot {snapshot_name}: {e}")