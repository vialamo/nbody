import numpy as np
from vispy import app, scene
from vispy.color import get_colormap

class SimulationApp:
    """Handles pure rendering and GUI events."""
    def __init__(self, engine):
        self.engine = engine  # Holds a reference to the backend physics engine
        self.config = engine.config
        
        self.canvas = scene.SceneCanvas(keys='interactive', 
                                        size=(self.config.render_size, self.config.render_size), 
                                        show=True, 
                                        title='N-Body + Hydrodynamics')
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = 'turntable' 
        self.view.camera.fov = 60
        self.view.camera.center = (self.config.domain_size/2, self.config.domain_size/2, self.config.domain_size/2)
        self.view.camera.distance = self.config.domain_size * 2

        try:
            self.colormap = get_colormap('plasma')
        except KeyError:
            # Fallback for older VisPy versions (like Ubuntu apt packages)
            print("Warning: 'plasma' colormap not found, falling back to 'hot'.")
            self.colormap = get_colormap('hot')
        
        if self.config.enable_hydro:
            initial_volume_data = np.zeros((self.config.mesh_size, self.config.mesh_size, self.config.mesh_size))
            self.gas_volume = scene.visuals.Volume(initial_volume_data, parent=self.view.scene, cmap=self.colormap, method='mip', threshold=0.1)
            self.gas_volume.transform = scene.transforms.STTransform(scale=(self.config.cell_size, self.config.cell_size, self.config.cell_size), translate=(0, 0, 0))
        
        particle_pos = np.zeros((len(self.engine.state.particles), 3), dtype=np.float32)
        self.particles_visual = scene.visuals.Markers(parent=self.view.scene)
        self.particles_visual.set_data(particle_pos, face_color='white', edge_color=None, size=2)
        
        box = scene.visuals.Box(width=self.config.domain_size, height=self.config.domain_size, depth=self.config.domain_size, color=(1, 1, 1, 0), edge_color=(1, 1, 1, 0.8), parent=self.view.scene)
        box.transform = scene.transforms.STTransform(translate=(self.config.domain_size/2, self.config.domain_size/2, self.config.domain_size/2))

        self.update_visuals() 

        # The VisPy timer now just tells the engine to run one step
        self.timer = app.Timer(interval='auto', connect=self.update_simulation, start=True)
        self.canvas.events.close.connect(self.on_close)

    def update_visuals(self):
        """Syncs the GPU visualization with the engine's current state."""
        if self.config.enable_hydro:
            log_field = np.log10(self.engine.state.gas.pressure + 1e-12)
            min_log, max_log = np.min(log_field), np.max(log_field)
            norm_field = (log_field - min_log) / (max_log - min_log) if max_log > min_log else np.zeros_like(log_field)
            volume_data = norm_field.transpose(2, 1, 0).astype(np.float32)
            self.gas_volume.set_data(volume_data)
            self.gas_volume.clim = (0.0, 1.0) 

        particle_pos = np.array([[p.x, p.y, p.z] for p in self.engine.state.particles], dtype=np.float32)
        self.particles_visual.set_data(particle_pos, face_color='white', edge_color='white', size=1)
        self.canvas.update()

    def update_simulation(self, event):
        """VisPy Event Callback."""
        # Step the physics forward
        self.engine.step()

        # Update visuals at a throttled rate (every 5 steps)
        if self.engine.cycle_count % 5 == 0:
            self.update_visuals()

    def on_close(self, event):
        self.timer.stop()
        self.canvas.app.quit()

    def run(self):
        self.canvas.app.run()