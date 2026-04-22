# Cosmological N-Body/Hydro Simulation Experiments

This repository documents my experiments in cosmological N-body/hydrodynamics simulations. It contains a "toy model" simulation program in C++, along with a book that explains the underlying physics and algorithms. This C++ implementation serves as a high-performance algorithmic testbed and uses Eigen and PocketFFT to explore memory-contiguous architectures and optimized linear algebra.

## Key Features Implemented

* **Gravity Solvers:**
    * **Particle-Particle (PP):** Direct summation for high-accuracy short-range forces.
    * **Particle-Mesh (PM):** FFT-based Poisson solver for efficient long-range forces.
    * **P³M (Particle-Particle Particle-Mesh):** A hybrid method combining PP and PM with a subtractive scheme and a smooth tapering function.
* **Cosmology:**
    * **Expanding Universe:** Simulation in comoving coordinates within an Einstein-de Sitter (EdS) model.
    * **Cosmological Integrator:** A Kick-Drift-Kick (KDK) Leapfrog scheme that correctly handles Hubble drag.
    * **Initial Conditions:** Particle generation on a lattice with perturbations applied via a simplified Zel'dovich Approximation using a power-law power spectrum.
    * **Adaptive Timestepping:** Dynamic calculation of the global timestep based on the Courant-Friedrichs-Lewy (CFL) hydro condition and maximum gravitational acceleration.
* **Hydrodynamics:**
    * **Grid-Based (Eulerian) Solver:** Implements a finite-volume solver for the adiabatic Euler equations on a fixed grid, tracking conservative variables (density, momentum, energy).
    * **HLL Riemann Solver:** Uses the Harten-Lax-van Leer (HLL) approximate Riemann solver to compute fluxes between cells.
    * **Operator Splitting:** Employs directional splitting (sequential X, Y, and Z-sweeps) to update the multidimensional grid.
    * **Two-Way Coupling:** The gas density contributes to the total gravitational field via the PM solver, and the gas momentum/energy is updated by gravitational source terms during the KDK kicks.
* **Numerical Methods:**
    * **Cloud-in-Cell (CIC):** A symmetric mass-assignment and force-interpolation scheme for the PM grid to ensure strict momentum conservation.
    * **Periodic Boundary Conditions:** A "wrap-around" universe to model a representative patch of a larger cosmos.
    * **Gravitational Softening:** Plummer softening to ensure numerical stability during close encounters.
* **Configuration & I/O:**
    * **Parameter Input:** All simulation parameters (domain size, particle count, physics, output timing, etc.) are read at runtime from a plain-text `simulation.ini` file.
    * **Snapshot Output:** Simulation snapshots (including particle data, gas grids, and metadata) are saved periodically using the HDF5 (Hierarchical Data Format) file format.
    * **Diagnostics:** Real-time logging of conservation metrics (mass, momentum, energy) and performance timings to stdout and CSV.
    
![N-Body Simulation Animation](simulation.gif)

## Repository Structure

* [`/src/`](src/): A high-performance C++ P³M + hydrodynamics cosmological simulation.
* [`/docs/`](docs/): Contains a book titled **"Notes on N-Body/Hydrodynamics Simulation"** in Markdown format (also in epub and pdf). This document summarizes the concepts, derivations, and algorithms implemented in the code.
* [`/scripts/`](scripts/): Contains Python utilities for post-processing, including a 3D visualizer (`viewer.py`) and a physical validation suite (`verify_run.py`) to analyze the HDF5 outputs.

## Getting Started

### C++ Cosmological code

1.  **Prerequisites (Linux/Ubuntu):**
    You need a C++ compiler, CMake, and the development libraries for HDF5, and Eigen.

    ```bash
    sudo apt update
    sudo apt install build-essential cmake libhdf5-dev libeigen3-dev libomp-dev
    ```

2.  **Compile (with CMake):**
    This project uses Modern CMake to find dependencies, download testing frameworks, and build the executables. The build is done in a separate `build` directory to keep the source folder clean.

    ```bash
    # From the project's root directory, configure the build system:
    cmake -B build
    cmake --build build
    ```

    This will create two executables inside the `build` directory: the main simulation `nbody` and the automated test suite `run_tests`.

3.  **Run:**
    The build process automatically copies the `simulation.ini` file into the `build` directory. You can run the simulation from there:

    ```bash
    # Navigate to the build directory:
    cd build
    # From inside the 'build' directory:
    ./nbody
    ```

    The program will start and automatically load its configuration from `simulation.ini` and create a uniquely timestamped folder inside an `outputs/` directory to safely store all snapshots and diagnostics.

4.  **Run the automated tests:**
    To ensure the physics engine and mathematical utilities are functioning correctly on your machine, you can run the built-in test suite:
    ```bash
    # From inside the 'build' directory:
    ctest --output-on-failure
    ```

### Viewer
1.  **Prerequisites:** Ensure you have Python 3 and the following libraries installed:
    ```bash
    pip install numpy vispy h5py PyOpenGL matplotlib pyqt5
    ```
    *(Note: VisPy requires a working OpenGL implementation on your system for the 3D volumetric rendering).*
2.  **Run:** Navigate to the `scripts/` directory (or wherever the script is). Run the script passing an output folder generated by a simulation as an argument:
    ```bash
    python viewer.py <folder>
    ```

### Physical Validation Tool
While the C++ `ctest` suite verifies the internal math during compilation, the `verify_run.py` script performs macroscopic physical validation on the generated HDF5 snapshots after a run. It checks for data integrity (NaNs, boundary escapes), strict thermodynamic positivity (preventing negative density/energy), and exact mass conservation across the entire simulation timeframe.

1.  **Prerequisites:** Requires `numpy` and `h5py` (which are already included if you installed the Viewer prerequisites).
2.  **Run:** Navigate to your project root and execute the script, passing the directory containing your `.hdf5` snapshots as an argument:
    ```bash
    python scripts/verify_run.py build/outputs/your_timestamped_folder/
    ```
    The script will print a report verifying the physical consistency of each snapshot.

## Learning Log & Guidebook

This project is developed as a learning exercise. The companion document, [`docs/book.pdf`](docs/book.pdf), is a "living book" that organizes and explains all the physics and computer science concepts encountered during this process. It is written in the style of a guide that I would have found most helpful when I began.
