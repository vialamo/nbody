import sys
import h5py
import numpy as np
from vispy import app, scene
from vispy.scene.visuals import Volume, Markers, Text
from vispy.color import get_colormap, Colormap
from vispy.io import write_png
import os

class True3DViewer:
    def __init__(self, filepath):
        self.filepath = filepath
        self.file = None
        
        # Initial load and configuration parsing
        self.refresh_file(initial=True)
        
        if self.num_frames == 0:
            print("Error: No snapshots found in the HDF5 file.")
            sys.exit(1)
            
        print(f"Loaded {self.filepath} - Found {self.num_frames} frames.")
        
        self.current_frame = 0
        self.is_playing = False
        
        self.setup_canvas()
        self.update_frame(0)

    def refresh_file(self, initial=False):
        """Handles opening the file, setting SWMR, and scanning for snapshots."""
        try:
            if self.file is not None:
                self.file.close()
            
            # Re-open safely with SWMR explicitly requested
            self.file = h5py.File(self.filepath, 'r', libver='latest', swmr=True)
            
            if initial:
                assert self.file.swmr_mode, "Failed to engage SWMR mode!"
                self.domain_size = self.file.attrs.get('domain_size', 1.0)
                self.use_hydro = bool(self.file.attrs.get('use_hydro', 0))
            
            # Re-read the available snapshots
            self.snapshots = sorted([k for k in self.file.keys() if k.startswith('snapshot_')])
            self.num_frames = len(self.snapshots)
            
            if not initial:
                print(f"Refreshed: Now seeing {self.num_frames} frames.")
                
        except Exception as e:
            if initial:
                raise e
            else:
                print(f"Warning: Refresh failed ({e}). Try again.")

    def setup_canvas(self):
        self.canvas = scene.SceneCanvas(keys='interactive', size=(1000, 800), show=True, bgcolor='black', title="Cosmological Render")
        self.view = self.canvas.central_widget.add_view()
        
        self.view.camera = scene.cameras.TurntableCamera(
            fov=45, 
            elevation=20, 
            azimuth=30, 
            distance=self.domain_size * 2.5
        )
        self.view.camera.center = (self.domain_size/2, self.domain_size/2, self.domain_size/2)

        # Create an Alpha-Gradient Colormap
        base_cmap = get_colormap('plasma')
        colors = base_cmap.map(np.linspace(0, 1, 256))
        colors[:, 3] = np.linspace(0, 1, 256) ** 2 
        self.custom_cmap = Colormap(colors)

        # Setup Gas Volume
        self.volume = None
        if self.use_hydro:
            self.volume = Volume(np.zeros((2, 2, 2)), cmap=self.custom_cmap, method='translucent', parent=self.view.scene)
            N = self.file.attrs.get('mesh_size', 64)
            scale_factor = self.domain_size / N
            self.volume.transform = scene.transforms.STTransform(scale=(scale_factor, scale_factor, scale_factor))

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
        print("[Right Arrow] : Next Frame")
        print("[Left Arrow]  : Previous Frame")
        print("[F5]          : Refresh File to load new frames")
        print("[E]           : Export all frames as clean PNGs")
        print("[Left Mouse]  : Rotate Camera")
        print("[Scroll]      : Zoom In/Out")
        print("[Shift + Left] : Pan Camera")

    def update_frame(self, frame_idx):
        group_name = self.snapshots[frame_idx]
        group = self.file[group_name]
        
        sim_time = group.attrs.get('simulation_time', 0.0)
        a = group.attrs.get('scale_factor', 1.0)
        z = (1.0 / a) - 1.0 if a > 0 else 0.0
        
        # Particles
        pos_x = group['particles/position_x'][:]
        pos_y = group['particles/position_y'][:]
        pos_z = group['particles/position_z'][:]
        positions = np.c_[pos_x, pos_y, pos_z]
        
        self.particles.set_data(positions, face_color=(1, 1, 1, 0.5), edge_width=0, size=2.0)

        # Gas
        p_max = 0
        if self.use_hydro and 'gas' in group:
            pressure = group['gas/pressure'][:]
            pressure = pressure.transpose(2, 1, 0)
            
            p_min, p_max = np.min(pressure), np.max(pressure)
            if p_max > p_min:
                normalized_pressure = (pressure - p_min) / (p_max - p_min)
                self.volume.set_data(normalized_pressure)
                self.volume.clim = [0.0, 1.0]

        # HUD
        hud_str = (
            f"Frame: {frame_idx:04d} / {self.num_frames-1}\n"
            f"Time:  {sim_time:.4f}\n"
            f"Scale: {a:.4f}\n"
            f"z:     {z:.2f}\n"
            f"Max P: {p_max:.2e}"
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
            self.refresh_file()

        elif event.key == 'E':
            self.export_frames()

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
        print("Usage: python viewer_vispy.py <path_to_hdf5_file>")
        sys.exit(1)
        
    viewer = True3DViewer(sys.argv[1])
    viewer.run()