import sys
import datetime
import h5py

from config import Config
from utils import Logger
from engine import SimulationEngine

def main():
    config_filename = sys.argv[1] if len(sys.argv) > 1 else "simulation.ini"
    config = Config(config_filename)
    
    run_timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = Logger(log_filename=f"sim_{run_timestamp}_diagnostics.csv")
    
    snapshot_filename = f"sim_{run_timestamp}.hdf5"
    h5_file = h5py.File(snapshot_filename, 'w')
    param_grp = h5_file.create_group("parameters")
    for key, val in config.get_all_params().items():
        param_grp.attrs[key] = val

    # Create the backend engine
    engine = SimulationEngine(config, logger, h5_file)

    try:
        if config.enable_visualization:
            from visualization import SimulationApp
            print("Starting graphical simulation loop...")
            sim_app = SimulationApp(engine)
            if sys.flags.interactive != 1:
                sim_app.run() # Script pauses here until window is closed
        else:
            print("Starting headless simulation loop...")
            while engine.cycle_count < config.max_cycles:
                engine.step()
                
    except KeyboardInterrupt:
        print("\nSimulation interrupted by user.")
        
    finally:
        print("Simulation finished. Closing file handles...")
        h5_file.close()
        logger.close()

if __name__ == "__main__":
    main()