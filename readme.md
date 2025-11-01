# Cosmological N-Body/Hydro Simulation Experiments

This repository documents my experiments in cosmological N-body/hydrodynamics simulations. It contains two "toy model" simulation programs, one in Python and one in C++, along with a book that explains the underlying physics and algorithms.

## Key Features Implemented

The simulations model a 2D universe and include a range of standard techniques used in cosmological codes:

* **Gravity Solvers:**
    * **Particle-Particle (PP):** Direct summation for high-accuracy short-range forces.
    * **Particle-Mesh (PM):** FFT-based Poisson solver for efficient long-range forces.
    * **P³M (Particle-Particle Particle-Mesh):** A hybrid method combining PP and PM with a subtractive scheme and a smooth tapering function.
* **Cosmology:**
    * **Expanding Universe:** Simulation in comoving coordinates within an Einstein-de Sitter (EdS) model.
    * **Cosmological Integrator:** A Kick-Drift-Kick (KDK) Leapfrog scheme that correctly handles Hubble drag.
    * **Initial Conditions:** Particle generation on a lattice with perturbations applied via a simplified Zel'dovich Approximation using a power-law power spectrum.
* **Hydrodynamics:**
    * **Grid-Based (Eulerian) Solver:** Implements a finite-volume solver for the adiabatic Euler equations on a fixed grid, tracking conservative variables (density, momentum, energy).
    * **HLL Riemann Solver:** Uses the Harten-Lax-van Leer (HLL) approximate Riemann solver to compute fluxes between cells.
    * **Operator Splitting:** Employs dimensional splitting (sequential X and Y-sweeps) to update the 2D grid.
    * **Two-Way Coupling:** The gas density contributes to the total gravitational field, and the gas itself is accelerated by gravity via source terms in the momentum and energy equations.
* **Numerical Methods:**
    * **Cloud-in-Cell (CIC):** A second-order mass assignment and force interpolation scheme for the PM grid.
    * **Periodic Boundary Conditions:** A "wrap-around" universe to model a representative patch of a larger cosmos.
    * **Gravitational Softening:** Plummer softening to ensure numerical stability during close encounters.
* **Configuration & I/O:**
    * **Parameter Input:** All simulation parameters (domain size, particle count, physics, output timing, etc.) are read at runtime from a plain-text `simulation.ini` file.
    * **Snapshot Output:** Simulation snapshots (including particle data, gas grids, and metadata) are saved periodically using the HDF5 (Hierarchical Data Format) file format.
    
![N-Body Simulation Animation](simulation.gif)

## Repository Structure

* [`/python/`](python/): A complete 2D P³M + hydrodynamics cosmological simulation written in Python.
* [`/cpp/`](cpp/): A high-performance C++ equivalent of the simulation.
* [`/book/`](book/): Contains a book titled **"Notes on N-Body/Hydrodynamics Simulation"** in Markdown format (also exported to epub and pdf). This document summarizes the concepts, derivations, and algorithms implemented in the code.

## Getting Started

### Python Version

1.  **Prerequisites:** Ensure you have Python 3 and the following libraries installed:
    ```bash
    pip install numpy pygame matplotlib scipy h5py
    ```
2.  **Run:** Navigate to the `python/` directory (or wherever the script is). Make sure the `simulation.ini` file is in the same directory. Run the script:
    ```bash
    python nbody.py
    ```

### C++ Version (Linux/Ubuntu)

1.  **Prerequisites:**
    You need a C++ compiler, CMake, and the development libraries for SFML, HDF5, and Eigen.

    ```bash
    sudo apt update
    sudo apt install build-essential cmake libsfml-dev libhdf5-dev libeigen3-dev
    ```

2.  **Compile (with CMake):**
    This project uses CMake to find all dependencies and build the executable. The build is done in a separate `build` directory to keep the source folder clean.

    ```bash
    # From the project's root directory (e.g., cpp/):
    mkdir build
    cd build

    # This finds all libraries and generates the Makefile
    cmake ..

    # This compiles the code into an executable
    make
    ```

    This will create an executable named `nbody` inside the `build` directory.

3.  **Run:**
    The build process automatically copies the `simulation.ini` file into the `build` directory. You can run the simulation from there:

    ```bash
    # From inside the 'build' directory:
    ./nbody
    ```

    The program will start and automatically load its configuration from `simulation.ini`.

## Learning Log & Guidebook

This project was developed as a learning exercise. The companion document, [`book/book.pdf`](book/book.pdf), is a "living book" that organizes and explains all the physics and computer science concepts encountered during this process. It is written in the style of a guide that I would have found most helpful when I began.
