import sys
import os
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt

import vispy
vispy.use('PyQt5') # Force the PyQt5 backend
from vispy import app, scene
from vispy.scene.visuals import Volume, Markers, Text
from vispy.color import Colormap
from vispy.io import write_png

def compute_cic_density(p_x, p_y, p_z, p_mass, domain_size, mesh_size):
    """Cloud-in-Cell mass splatting to estimate density."""
    N = mesh_size
    grid = np.zeros((N, N, N), dtype=np.float32)
    
    x = (p_x / domain_size) * N
    y = (p_y / domain_size) * N
    z = (p_z / domain_size) * N
    
    ix, iy, iz = np.floor(x).astype(int) % N, np.floor(y).astype(int) % N, np.floor(z).astype(int) % N
    fx, fy, fz = x - np.floor(x), y - np.floor(y), z - np.floor(z)
    ix1, iy1, iz1 = (ix + 1) % N, (iy + 1) % N, (iz + 1) % N
    
    np.add.at(grid, (ix, iy, iz), (1-fx)*(1-fy)*(1-fz) * p_mass)
    np.add.at(grid, (ix1, iy, iz), fx*(1-fy)*(1-fz) * p_mass)
    np.add.at(grid, (ix, iy1, iz), (1-fx)*fy*(1-fz) * p_mass)
    np.add.at(grid, (ix1, iy1, iz), fx*fy*(1-fz) * p_mass)
    
    np.add.at(grid, (ix, iy, iz1), (1-fx)*(1-fy)*fz * p_mass)
    np.add.at(grid, (ix1, iy, iz1), fx*(1-fy)*fz * p_mass)
    np.add.at(grid, (ix, iy1, iz1), (1-fx)*fy*fz * p_mass)
    np.add.at(grid, (ix1, iy1, iz1), fx*fy*fz * p_mass)
    
    return grid

class True3DViewer:
    def __init__(self, run_dir):
        self.run_dir = run_dir
        
        # Initial load and configuration parsing
        self.refresh_directory(initial=True)
        
        if self.num_frames == 0:
            print(f"Error: No snapshot files found in '{self.run_dir}'.")
            sys.exit(1)
            
        print(f"Loaded directory '{self.run_dir}' - Found {self.num_frames} frames.")
        
        self.current_frame = 0
        self.is_playing = False
        self.fixed_colormap = True
        
        self.setup_canvas()
        self.update_frame(0)

    def refresh_directory(self, initial=False):
        """Scans the directory for snapshot files."""
        try:
            search_pattern = os.path.join(self.run_dir, "snapshot_*.hdf5")
            self.snapshot_files = sorted(glob.glob(search_pattern))
            self.num_frames = len(self.snapshot_files)
            
            if initial and self.num_frames > 0:
                # Extract static simulation parameters from the first snapshot
                with h5py.File(self.snapshot_files[0], 'r', libver='latest', swmr=True) as f:
                    config_attr = f['Config'].attrs
                    self.domain_size = config_attr.get('domain_size', 1.0)
                    self.use_hydro = bool(config_attr.get('use_hydro', 0))
                    mesh_size = config_attr.get('mesh_size', 64)
                
                print("Reading the final snapshot to establish stable colormap bounds...")
                with h5py.File(self.snapshot_files[-1], 'r', libver='latest', swmr=True) as f:
                    # Pre-calculate Gas Extents
                    if self.use_hydro and 'Gas' in f:
                        if 'Gas/temperature' in f:
                            temp = f['Gas/temperature'][:]
                        else:
                            pressure = f['Gas/pressure'][:]
                            density = f['Gas/density'][:]
                            KB, M_H = 1.380649e-16, 1.6726219e-24
                            MU = f['Config'].attrs.get('primordial_mu', 1.22)
                            v_unit = f['Units'].attrs.get('unit_velocity_in_cgs', 1.0)
                            temp = (pressure / density) * (v_unit**2) * (MU * M_H / KB)
                            
                        self.global_t_max_log = np.max(np.log10(np.clip(temp, 10.0, None)))
                        self.global_t_min_log = 1.0  # Safe floor of 10 K
                    else:
                        self.global_t_min_log, self.global_t_max_log = 1.0, 5.0

                    # Pre-calculate DM Density Extents
                    if 'Particles/position_x' in f:
                        pos_x = f['Particles/position_x'][:]
                        pos_y = f['Particles/position_y'][:]
                        pos_z = f['Particles/position_z'][:]
                        mass  = f['Particles/mass'][:]
                        
                        dm_grid = compute_cic_density(pos_x, pos_y, pos_z, mass, self.domain_size, mesh_size)
                        safe_dens = np.clip(dm_grid, 1e-10, None)
                        self.global_dm_max_log = np.max(np.log10(safe_dens))
                        self.global_dm_min_log = np.min(np.log10(safe_dens))
                    else:
                        self.global_dm_min_log, self.global_dm_max_log = -10.0, 0.0
            
            if not initial:
                print(f"Refreshed: Now seeing {self.num_frames} frames.")
                
        except Exception as e:
            if initial:
                raise e
            else:
                print(f"Warning: Refresh failed ({e}). Try again.")

    def update_colorbar_labels(self, t_min_log, t_max_log):
        """Update the colorbar labels with current temperature range."""
        if hasattr(self, 'cbar_min'):
            self.cbar_min.text = f"10e{t_min_log:.1f} K"
            self.cbar_max.text = f"10e{t_max_log:.1f} K"

    def setup_with_colorbar(self):
        """Setup canvas with a separate colorbar widget."""
        from vispy.scene.widgets import Grid, ViewBox
        
        # Create a grid layout
        self.grid = self.canvas.central_widget.add_grid()
        
        # Main 3D view (takes most space)
        self.view = self.grid.add_view(row=0, col=1, col_span=10)
        self.view.camera = scene.cameras.TurntableCamera(
            fov=45, elevation=20, azimuth=30, distance=self.domain_size * 2.5
        )
        self.view.camera.center = (self.domain_size/2, self.domain_size/2, self.domain_size/2)
        
        # Colorbar view
        cbar_view = self.grid.add_view(row=0, col=0, col_span=1)
        cbar_view.camera = scene.cameras.PanZoomCamera(aspect=0.3, rect=(0, 0, 1, 1))
        
        # Create the colorbar visual
        n_colors = 256
        colors = plt.get_cmap('plasma')(np.linspace(0, 1, n_colors))
        
        # Use a series of rectangles
        rect_height = 1.0 / n_colors
        
        for i in range(n_colors):
            rect = scene.visuals.Rectangle(
                center=(0.5, i * rect_height + rect_height/2),
                width=0.2,
                height=rect_height,
                color=colors[i],
                parent=cbar_view.scene
            )
        
        # Add labels
        self.cbar_min = Text("10³ K", color='white', pos=(0.5, 0.0), anchor_x='center', 
                             anchor_y='bottom', parent=cbar_view.scene, font_size=8)
        self.cbar_max = Text("10⁷ K", color='white', pos=(0.5, 1.0), anchor_x='center', 
                             anchor_y='top', parent=cbar_view.scene, font_size=8)

    def setup_canvas(self):
        self.canvas = scene.SceneCanvas(keys='interactive', size=(1000, 800), show=True, 
                                        bgcolor='black', title="Cosmological Simulation Viewer")

        self.setup_with_colorbar()

        # Create an Alpha-Gradient Colormap
        colors = plt.get_cmap('plasma')(np.linspace(0, 1, 256))
        
        # Force the Alpha (opacity) channel to fade quadratically
        colors[:, 3] = np.linspace(0, 1, 256) ** 2
        colors[0] = [0.0, 0.0, 0.0, 0.0]
        self.custom_cmap = Colormap(colors)

        # Setup Gas Volume
        self.volume = None
        if self.use_hydro:
            self.volume = Volume(np.zeros((2, 2, 2), dtype=np.float32), cmap=self.custom_cmap, 
                                 method='translucent', parent=self.view.scene)
            # Read mesh size from the first file to scale the volume properly
            with h5py.File(self.snapshot_files[0], 'r') as f:
                N = f['Config'].attrs.get('mesh_size', 64)
            scale_factor = self.domain_size / N
            self.volume.transform = scene.transforms.STTransform(
                scale=(scale_factor, scale_factor, scale_factor))

        # Setup Dark Matter Particles
        self.particles = Markers(parent=self.view.scene)
        self.particles.set_gl_state('translucent', blend=True, depth_test=False)

        # Setup HUD
        self.hud_text = Text(
            text='', color='white', pos=(20, 20), anchor_x='left', anchor_y='bottom',
            font_size=12, parent=self.canvas.scene 
        )

        self.holding_right = False
        self.holding_left = False

        # Playback Timer & Events
        self.timer = app.Timer('auto', connect=self.on_timer, start=False)
        self.canvas.events.key_press.connect(self.on_key_press)
        self.canvas.events.key_release.connect(self.on_key_release)
        
        print("\n--- Controls ---")
        print("[Spacebar] : Play / Pause")
        print("[Right Arrow] : Next snapshot")
        print("[Left Arrow]  : Previous snapshot")
        print("[F5]          : Refresh folder to load new snapshots")
        print("[E]           : Export all frames as clean PNGs")
        print("[C]           : Toggle fixed / dynamic colormap")
        print("[Left Mouse]  : Rotate camera")
        print("[Scroll]      : Zoom in/out")
        print("[Shift + Left] : Pan camera")

    def update_frame(self, frame_idx):
        filepath = self.snapshot_files[frame_idx]
        
        with h5py.File(filepath, 'r', libver='latest', swmr=True) as f:
            header_attr = f['Header'].attrs
            config_attr = f['Config'].attrs
            units_attr  = f['Units'].attrs

            # Pull dynamic telemetry from /Header
            sim_time = header_attr.get('simulation_time', 0.0)
            a = header_attr.get('scale_factor', 1.0)
            z = (1.0 / a) - 1.0 if a > 0 else 0.0

            unit_time_gyr = units_attr.get('unit_time_in_gyr', 1.0)
            unit_length_mpc = units_attr.get('unit_length_in_mpc', 1.0)
            
            # Calculate physical time and sizes
            time_gyr = sim_time * unit_time_gyr
            comoving_box_mpc = self.domain_size * unit_length_mpc
            physical_box_mpc = comoving_box_mpc * a           

            # Particles
            pos_x = f['Particles/position_x'][:]
            pos_y = f['Particles/position_y'][:]
            pos_z = f['Particles/position_z'][:]
            mass  = f['Particles/mass'][:]
            positions = np.c_[pos_x, pos_y, pos_z]
            
            # Estimate local density using the CIC grid
            N = config_attr.get('mesh_size', 64)
            dm_grid = compute_cic_density(pos_x, pos_y, pos_z, mass, self.domain_size, N)
            
            # Find which grid cell each particle is currently in
            ix = np.floor((pos_x / self.domain_size) * N).astype(int) % N
            iy = np.floor((pos_y / self.domain_size) * N).astype(int) % N
            iz = np.floor((pos_z / self.domain_size) * N).astype(int) % N
            
            # Sample the local density for each particle (Nearest Grid Point)
            particle_densities = dm_grid[ix, iy, iz]
            
            # Convert to Logarithmic Opacity
            safe_dens = np.clip(particle_densities, 1e-10, None)
            log_dens = np.log10(safe_dens)
            
            if self.fixed_colormap:
                min_log = self.global_dm_min_log
                max_log = self.global_dm_max_log
            else:
                min_log = np.min(log_dens)
                max_log = np.max(log_dens)
            
            # Create the RGBA color array
            colors = np.zeros((len(pos_x), 4), dtype=np.float32)
            colors[:, 0] = 1  # R
            colors[:, 1] = 1  # G
            colors[:, 2] = 1  # B
            
            if max_log > min_log:
                alphas = (log_dens - min_log) / (max_log - min_log)
                # Apply a gamma curve to hide low-density void particles
                gamma_curve = 3.0 
                alphas = alphas ** gamma_curve
                
                # Cap max opacity so halo cores still look translucent
                colors[:, 3] = alphas * 0.65 
            else:
                colors[:, 3] = 0.1

            # Update the markers with the custom per-particle color array
            self.particles.set_data(positions, face_color=colors, edge_width=0, size=1.25)

            # Gas
            t_max = 0
            if self.use_hydro and 'Gas' in f:
                if 'Gas/temperature' in f:
                    temperature = f['Gas/temperature'][:]
                    temperature = temperature.transpose(2, 1, 0)
                else:
                    # Fallback for adiabatic runs without the temperature
                    pressure = f['Gas/pressure'][:].transpose(2, 1, 0)
                    density = f['Gas/density'][:].transpose(2, 1, 0)
                    
                    KB = 1.380649e-16
                    M_H = 1.6726219e-24
                    MU = config_attr.get('primordial_mu', 1.22)
                    v_unit_cgs = units_attr.get('unit_velocity_in_cgs', 1.0)
                    
                    specific_energy_cgs = (pressure / density) * (v_unit_cgs**2)
                    temperature = specific_energy_cgs * (MU * M_H / KB)

                t_max = np.max(temperature)
                
                # Prevent log(0) warnings by clipping to the 10.0K floor
                safe_temp = np.clip(temperature, 10.0, None) 
                log_temp = np.log10(safe_temp)
                            
                if self.fixed_colormap:
                    t_min_log = self.global_t_min_log
                    t_max_log = self.global_t_max_log
                else:
                    t_min_log, t_max_log = np.min(log_temp), np.max(log_temp)
                
                if t_max_log > t_min_log:
                    normalized_temp = (log_temp - t_min_log) / (t_max_log - t_min_log)
                    gamma = 2.5
                    normalized_temp = normalized_temp ** gamma
                    normalized_temp = normalized_temp.astype(np.float32)
                    self.volume.set_data(normalized_temp)
                    self.volume.clim = [0.0, 1.0]

        if (not hasattr(self, 'temp_min_log') 
            or self.temp_min_log != t_min_log 
            or self.temp_max_log != t_max_log):
            self.temp_min_log = t_min_log
            self.temp_max_log = t_max_log
            self.update_colorbar_labels(t_min_log, t_max_log)

        # HUD
        hud_str = (
            f"Snapshot: {frame_idx:04d} / {self.num_frames-1}\n"
            f"Epoch: a = {a:.4f}  |  z = {z:.2f}\n"
            f"Time:  {sim_time:.4f} code  ({time_gyr:.2f} Gyr)\n"
            f"Size:  {comoving_box_mpc:.1f} Mpc comoving  ({physical_box_mpc:.1f} Mpc physical)\n"
            f"Max T: {t_max:.2e} K"
        )
        self.hud_text.text = hud_str
        self.canvas.update()

    def export_frames(self):
        """Exports all frames as clean PNGs without the HUD."""
        export_dir = "export_frames"
        os.makedirs(export_dir, exist_ok=True)
        
        print(f"\nStarting export of {self.num_frames} frames to '{export_dir}/'...")
        
        # Save current state (user's session)
        original_frame = self.current_frame
        was_playing = self.is_playing
        if self.is_playing:
            self.timer.stop()
            self.is_playing = False
            
        for i in range(self.num_frames):
            self.update_frame(i)
            self.hud_text.text = ''
            app.process_events() 
            
            img = self.canvas.render()
            filepath = os.path.join(export_dir, f"frame_{i:04d}.png")
            write_png(filepath, img)
            
            sys.stdout.write(f"\rExporting: {i+1}/{self.num_frames} ({filepath})")
            sys.stdout.flush()
            
        print("\nExport complete!")
        
        # Restore the viewer back
        self.update_frame(original_frame)
        if was_playing:
            self.timer.start()
            self.is_playing = True

    def on_timer(self, event):
        if self.holding_right:
            if self.current_frame < self.num_frames - 1:
                self.current_frame += 1
                self.update_frame(self.current_frame)
                
        elif self.holding_left:
            if self.current_frame > 0:
                self.current_frame -= 1
                self.update_frame(self.current_frame)
                
        elif self.is_playing:
            if self.current_frame < self.num_frames - 1:
                self.current_frame += 1
                self.update_frame(self.current_frame)
            else:
                self.timer.stop()
                self.is_playing = False

    def on_key_press(self, event):
        if event.key == 'Space':
            self.is_playing = not self.is_playing
            if self.is_playing:
                if self.current_frame >= self.num_frames - 1:
                    self.current_frame = 0
                self.timer.start()
            else:
                if not self.holding_left and not self.holding_right:
                    self.timer.stop()
                    
        elif event.key == 'Right':
            self.holding_right = True
            self.timer.start()
            
        elif event.key == 'Left':
            self.holding_left = True
            self.timer.start()
            
        elif event.key == 'F5':
            self.refresh_directory()

        elif event.key == 'E':
            self.export_frames()
        elif event.key == 'C':
            self.fixed_colormap = not self.fixed_colormap
            self.update_frame(self.current_frame) # Force refresh

    def on_key_release(self, event):
        if event.key == 'Right':
            self.holding_right = False
            if not self.is_playing and not self.holding_left:
                self.timer.stop()
                
        elif event.key == 'Left':
            self.holding_left = False
            if not self.is_playing and not self.holding_right:
                self.timer.stop()

    def run(self):
        app.run()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python viewer_vispy.py <path_to_run_directory>")
        sys.exit(1)
        
    viewer = True3DViewer(sys.argv[1])
    viewer.run()