# Cosmological Simulations

This repository documents my experiments in cosmological N-body/hydrodynamics simulations. It contains a simulation program in C++, along with a book that explains the underlying physics and algorithms. This C++ implementation serves as a high-performance algorithmic testbed to explore memory-contiguous architectures and optimized linear algebra.

## Key Features Implemented

* **Gravity Solvers:**
    * **Particle-Particle (PP):** Direct summation for high-accuracy short-range forces.
    * **Particle-Mesh (PM):** FFT-based Poisson solver for efficient long-range forces.
    * **P³M (Particle-Particle Particle-Mesh):** A hybrid method combining PP and PM with a subtractive scheme and a smooth tapering function.
* **Cosmology:**
    * **Expanding Universe:** Simulation in comoving coordinates within an Einstein-de Sitter (EdS) model.
    * **Cosmological Integrator:** A Kick-Drift-Kick (KDK) Leapfrog scheme that correctly handles Hubble drag.
    * **Initial Conditions:** Advanced cosmological initial conditions. Particles are perturbed from a uniform lattice using the full 3D Zel'dovich Approximation, deriving physical displacements and velocities from a Gaussian random field generated in Fourier space via a cosmological power spectrum.
    * **Adaptive Timestepping:** Dynamic calculation of the global timestep based on the Courant-Friedrichs-Lewy (CFL) hydro condition and maximum gravitational acceleration.
* **Hydrodynamics:**
    * **Grid-Based (Eulerian) Solver:** Implements a finite-volume solver for the adiabatic Euler equations on a fixed grid, tracking conservative variables (density, momentum, energy).
    * **HLL Riemann Solver:** Uses the Harten-Lax-van Leer (HLL) approximate Riemann solver to compute fluxes between cells.
    * **Operator Splitting:** Employs directional splitting (sequential X, Y, and Z-sweeps) to update the multidimensional grid.
    * **Two-Way Coupling:** The gas density contributes to the total gravitational field via the PM solver, and the gas momentum/energy is updated by gravitational source terms during the KDK kicks.
* **High-Performance Computing (HPC):**
    * **OpenMP Multithreading:** Heavy loops (such as the Riemann solver, mass assignment, and grid calculations) are parallelized across available CPU cores.
    * **SIMD Vectorization:** Core CPU mathematics leverage Eigen and AVX vectorization for cache-friendly, contiguous memory speedups.
    * **GPU Offloading:** The computationally expensive direct-summation particle forces can be optionally offloaded to NVIDIA GPUs using OpenMP `#pragma omp target` directives and the Clang/LLVM toolchain.
* **Numerical Methods:**
    * **Cloud-in-Cell (CIC):** A symmetric mass-assignment and force-interpolation scheme for the PM grid to ensure strict momentum conservation.
    * **Periodic Boundary Conditions:** A "wrap-around" universe to model a representative patch of a larger cosmos.
    * **Gravitational Softening:** Plummer softening to ensure numerical stability during close encounters.
* **Configuration & I/O:**
    * **Parameter Input:** All simulation parameters (domain size, particle count, physics, output timing, etc.) are read at runtime from a plain-text `simulation.ini` file.
    * **Snapshot Output:** Simulation snapshots (including particle data, gas grids, and metadata) are saved periodically using the HDF5 (Hierarchical Data Format) file format.
    * **Diagnostics:** Real-time logging of conservation metrics (mass, momentum, energy) and performance timings to stdout and CSV.
    
![Cosmological Simulation Animation](simulation.gif)

## Repository Structure

* [`/src/`](src/): A high-performance C++ P³M + hydrodynamics cosmological simulation.
* [`/docs/`](docs/): Contains a "living book" titled **"Notes on Cosmological Simulations"** in Markdown format. It includes a `Makefile` to automatically build the source notes into EPUB and PDF files.
* [`/scripts/`](scripts/): Contains Python utilities for post-processing, including a 3D visualizer (`viewer.py`) and a physical validation suite (`verify_run.py`) to analyze the HDF5 outputs.

## Getting Started

### C++ Cosmological code
This section refers to the ASIMOV (Advanced Simulation of Intergalactic Matter and Observable Voids) code.

1.  **Prerequisites (Linux/Ubuntu):**
    You need a C++ compiler, CMake, and the development libraries for HDF5, Eigen, and OpenMP. To utilize GPU offloading, the Clang/LLVM toolchain is also required.

    ```bash
    sudo apt update
    sudo apt install build-essential cmake libhdf5-dev libeigen3-dev libomp-dev
    # Optional: Install Clang for GPU Offloading support
    sudo apt install clang lld
    ```

2.  **Compile (with CMake):**
    This project uses Modern CMake to find dependencies, download testing frameworks, and build the executables. You can build the project for standard CPU execution (using AVX/SIMD) or enable NVIDIA GPU offloading.

    **Option A: Standard CPU Build (Default)**
    ```bash
    mkdir build && cd build
    cmake ..
    make -j
    ```

    **Option B: GPU Offload Build (Requires Clang & NVIDIA GPU)**
    ```bash
    mkdir build && cd build
    # Force CMake to use Clang and enable the GPU offload toggle
    CXX=clang++ CC=clang cmake -DUSE_GPU_OFFLOAD=ON ..
    make -j
    ```

    This will create two executables inside the `build` directory: the main simulation `asimov` and the automated test suite `run_tests`.

3.  **Run:**
    The build process automatically copies the `simulation.ini` file into the `build` directory. You can run the simulation from there:

    ```bash
    # From inside the 'build' directory:
    ./asimov
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
    python scripts/verify_run.py build/outputs/timestamped_folder/
    ```
    The script will print a report verifying the physical consistency of each snapshot.

### Building the Book (Documentation)
The `/docs` folder contains a guide detailing the physics and algorithms used in this project. You can compile the Markdown source into EPUB and PDF documents.

1.  **Prerequisites:** You need Pandoc and a LaTeX distribution to generate the PDF.
    ```bash
    sudo apt update
    sudo apt install pandoc texlive-latex-base texlive-latex-extra
    ```
2.  **Build:** Navigate to the `docs` directory and use the provided `Makefile`.
    ```bash
    cd docs
    make        # Builds both EPUB and PDF
    
    # Or build them individually:
    make pdf
    make epub
    ```

## Learning Log & Guidebook

This project is developed as a learning exercise. The companion document, built from the source files in [`/docs/`](docs/), is a "living book" that organizes and explains all the physics and computer science concepts encountered during this process. It is written in the style of a guide that I would have found most helpful when I began.

## Configuration Parameters

The simulation is configured using a standard INI file format. The configuration dictates the physical properties of the simulated universe, the grid resolution, the integration timesteps, and the data output frequency. Below is a breakdown of all available parameters by section.

### `[domain]`

Defines the spatial properties and resolution of the simulation box.

* **`domain_size`**: The internal code-unit length of the box (typically set to `1.0` for natural units).
* **`mesh_size`**: The number of grid cells along one axis for the Poisson solver and hydrodynamics (e.g., `32` creates a 32x32x32 computational grid).
* **`box_size_mpc`**: The physical comoving size of the simulation box in Megaparsecs (Mpc). Used to map internal code units back to physical reality.

### `[cosmology]`

Defines the cosmological model, specifically the energy budget and expansion rate of the universe.

* **`omega_baryon`**: The density parameter for normal (baryonic) matter.
* **`omega_M`**: The total matter density parameter (baryons + dark matter). Dark energy is dynamically assumed to be `1.0 - omega_M` for a flat universe.
* **`omega_lambda`**: The dark energy density parameter (cosmological constant).
* **`hubble_param`**: The dimensionless physical Hubble parameter (*h*), where H0 = 100 * h km/s/Mpc.
* **`expanding_universe`**: Boolean (`true` or `false`). If true, the code applies Hubble drag and stretches the comoving grid over time.

### `[initial_conditions]`

Controls the generation of the primordial density field and particle distributions.

* **`spectral_index`**: The primordial power spectrum index (*n_s*), typically `0.96` based on Planck data.
* **`start_a`**: The scale factor at which the simulation begins (e.g., `0.02` corresponds to redshift z = 49).
* **`sigma_8`**: The normalization of the power spectrum, defining the amplitude of density fluctuations at an 8 Mpc/h scale.
* **`initial_gas_temp_k`**: The physical temperature of the baryonic gas at `start_a`, in Kelvin.
* **`n_per_side`**: The number of N-body particles along one axis (e.g., `32` creates 32,768 total particles). Often matched to `mesh_size` to avoid interpolation artifacts.
* **`standing_particles`**: Boolean. If `true`, all N-body particles remain completely stationary throughout the entire duration of the simulation. Their positions and velocities are never updated. This is strictly a debugging feature, useful for isolating and testing other components of the code, such as observing how the hydrodynamics solver behaves around a static gravitational potential.
* **`seed`**: Integer seed for the random number generator, ensuring reproducible initial density fields.

### `[hydro]`

Configures the fluid dynamics solver for the baryonic gas.

* **`use_hydro`**: Boolean flag to enable or disable the hydrodynamics solver. If `false`, the code runs as a dark-matter-only N-body simulation.
* **`gamma`**: The adiabatic index (ratio of specific heats) of the gas. Set to `1.6666666667` (5/3) for a monatomic, non-relativistic ideal gas.

### `[p3m]`

Configures the Particle-Particle Particle-Mesh (P³M) gravity solver.

* **`use_pm`**: Boolean. Enables the long-range Particle-Mesh (PM) force calculation via Fast Fourier Transform.
* **`use_pp`**: Boolean. Enables the short-range Particle-Particle (PP) direct summation for sub-grid resolution.
* **`cutoff_radius_cells`**: The matching radius (*r_c*) where the algorithm switches from short-range PP forces to long-range PM forces, defined in units of grid cells.
* **`cutoff_transition_width_factor`**: Determines the width of the smoothing kernel used to seamlessly blend the PP and PM forces at the cutoff boundary to avoid unphysical jumps in acceleration.

### `[time]`

Controls the adaptive timestepping and integration limits.

* **`dt_factor`**: A global multiplier for the gravitational timestep, acting as a fraction of the system's local dynamical time.
* **`cfl_safety_factor`**: The Courant-Friedrichs-Lewy (CFL) number. Restricts the hydrodynamics timestep to ensure information does not travel further than one grid cell per step (must be < 1.0, typically `0.4`).
* **`max_scale_factor`**: The scale factor at which the simulation terminates. Set to `1.0` to run up to the present day (z = 0).

### `[output]`

Manages how and when the simulation writes data to disk.

* **`save_hdf5_every_delta_a`**: The interval for writing full snapshot files (particles, mesh densities, velocities) to disk, measured in scale factor increments.
* **`debug_info_every_cycles`**: The frequency (in integration steps) at which the code prints its current status, timestep, and performance metrics to the console/log.
* **`enable_energy_diagnostics`**: Boolean. If `true`, the code continuously calculates and verifies the conservation of energy and momentum, writing the error margins to a diagnostic file.