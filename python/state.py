class SimState:
    """Holds the entire physical state of the simulation."""
    def __init__(self, particles, gas, total_time=0.0, scale_factor=1.0):
        self.particles = particles
        self.gas = gas
        self.total_time = total_time
        self.scale_factor = scale_factor