# A Hybrid PM+PP/MFM Cosmological Hydrodynamics Code with GPU Acceleration



## Overview

We present a new cosmological hydrodynamics code combining a hybrid Particle-Mesh + Particle-Particle (PM+PP) gravity solver with a Meshless Finite Mass (MFM) hydrodynamics scheme. The code is written in C++17 with OpenMP thread parallelism and optional GPU acceleration via OpenMP target offload, targeting modern heterogeneous HPC architectures. Gravity is computed using a split PM+PP approach with a Fourier-space Poisson solver for the long-range force and an exact short-range PP correction accelerated by a Linear Bounding Volume Hierarchy (LBVH) tree constructed using the Karras (2012) algorithm with Morton-encoded 63-bit keys. Hydrodynamics follows the MFM formulation of Hopkins (2015), using a cubic-spline kernel, least-squares gradient estimation with Barth-Jespersen and pairwise slope limiting, and an HLLC Riemann solver with an exact Newton-Raphson runtime alternative. The code supports flat $\Lambda$CDM and Einstein-de Sitter cosmologies, individual particle timesteps with the Saitoh & Makino (2009) timestep limiter, a dual-energy formalism with an entropy switch, and radiative cooling via Grackle-format Haardt & Madau (2012) UV-background tables. Initial conditions can be generated internally using the Zel'dovich (1LPT) approximation with a BBKS transfer function normalized to $\sigma_8$. We describe the code design, numerical methods, and validate against standard hydrodynamic benchmarks including the Sod shock tube, Sedov–Taylor blast wave, linear sound wave, Kelvin–Helmholtz instability, Evrard collapse, and cosmological Zel'dovich pancake test.



## 1. Introduction

Production cosmological gravity solvers generally fall into three categories: adaptive mesh refinement (AMR) multigrid solvers (e.g., RAMSES, Teyssier 2002; Enzo, Bryan et al. 2014), pure hierarchical methods such as the Fast Multipole Method (FMM; e.g., pkdgrav3, Potter, Stadel & Teyssier 2017), and hybrid Particle-Mesh (PM) schemes that pair an FFT-based long-range solver with a short-range correction. Hybrid methods are particularly widespread—used in Gadget-2 (Springel 2005), Gadget-4 (Springel et al. 2021), and SWIFT (Schaller et al. 2024) via TreePM or FMM-PM splits—because they naturally handle periodic boundary conditions while isolating short-range forces to a compact cutoff radius. Adapting these short-range solvers to modern GPU architectures, however, remains challenging, as pointer-based octrees and high-order multipole expansions introduce thread divergence.

For hydrodynamics, the community has converged on three broad paradigms: Smoothed Particle Hydrodynamics (SPH, Monaghan 1992; Price 2012), adaptive mesh refinement (AMR, Teyssier 2002; Bryan et al. 2014), and meshless finite-volume/mass methods. AREPO (Springel 2010) introduced the moving-mesh approach; GIZMO (Hopkins 2015) introduced Meshless Finite Mass (MFM) and Meshless Finite Volume (MFV) methods that combine advantages of both particle-based and grid-based approaches, capturing shocks accurately without the artificial viscosity inherent to standard SPH or the mesh-motion overhead of AREPO.

Here we present a new code developed as a testbed for evaluating alternative, GPU-friendly algorithms on heterogeneous hardware. Our approach uses a PM+PP gravity split with an LBVH tree—a flat, cache-coherent tree structure originally developed for ray-tracing workloads (Lauterbach et al. 2009; Karras 2012)—combined with the MFM hydrodynamics scheme of Hopkins (2015). The code supports optional GPU execution of the LBVH construction and short-range PP gravity kernel via OpenMP target offloading, enabling portability across vendor GPU backends, as well as concurrent CPU hydrodynamics and GPU gravity execution.

This document is organized as follows. Section 2 describes the overall code architecture. Section 3 covers the cosmological framework and unit system. Section 4 describes initial conditions generation. Sections 5 and 6 detail the gravity and LBVH tree solvers. Section 7 covers MFM hydrodynamics. Section 8 addresses thermodynamics and the dual-energy formalism. Section 9 covers radiative cooling. Section 10 describes the time integration scheme. Section 11 presents code validation tests. Section 12 discusses performance. Section 13 covers data output. We discuss future work in Section 14.



## 2. Code Architecture

### 2.1 Language and Dependencies

The code is written in C++17 and relies on the following libraries:

- **OpenMP 5.0+**: CPU thread parallelism and GPU offload.
- **FFTW3** (Frigo & Johnson 2005): 3D real-to-complex FFTs for the Particle-Mesh solver and cosmological initial conditions, with multi-threaded support.
- **Eigen 3** (Guennebaud et al. 2010): dense linear algebra for matrix inversion, eigendecomposition, and vector operations used in gradient estimation and reconstruction.
- **HDF5 / H5Cpp** (The HDF Group): snapshot I/O in a custom HDF5 schema with Single-Writer Multiple-Reader (SWMR) mode for live inspection during runs, as well as reading of glass initial condition files and cooling tables.

### 2.2 Memory Layout

Particle data are stored in Structure-of-Arrays (SoA) layout, with dark matter and gas held in separate containers. This layout favors SIMD vectorization of inner loops and GPU memory coalescing. Gas particles carry additional fields: smoothing length $h$, density $\rho$, pressure $P$, specific internal energy $u$, specific total energy $E$, entropic function $A$, metallicity fraction, least-squares inverse matrices $\mathbf{B}_i$, primitive gradients ($\nabla\rho$, $\nabla P$, $\nabla v_x$, $\nabla v_y$, $\nabla v_z$), hydrodynamic impulse accumulators, predicted primitive states, and dual-energy/timestep tracking variables.

### 2.3 Parallelism Model

Shared-memory parallelism is obtained via OpenMP, using dynamic scheduling for irregular neighbor-search kernels (density, hydrodynamics, CPU PP gravity) and static scheduling for regular loops (integrator, PM, cooling). When compiled with GPU support, the PP gravity traversal offloads to the GPU via OpenMP target directives. Unified DM+gas arrays and the LBVH are kept resident in device memory, executing full GPU LBVH construction and traversal. GPU gravity and CPU hydrodynamics are run in parallel by default. When compiled without GPU support, the code falls back to the CPU path.

### 2.4 Control Flow

The main loop has the following structure:

```
while (a < a_end) {
    1. Sync clocks & process Saitoh-Makino wake-ups
    2. Kick 1          (half-step gravity + Hubble drag + thermodynamics)
    3. Drift           (full-step comoving position update)
    4. Advance clock & update cosmology (a, H)
    5. Update softening & rebuild LBVH tree(s)
    6. Compute long-range PM gravity
    7. Compute gas density, smoothing lengths h, and B matrices
    8. Compute MFM hydro (predict primitives, gradients, Riemann fluxes)
       & short-range PP gravity (optionally overlapped on CPU/GPU)
    9. Kick 2          (half-step gravity + hydro impulses + DEF/entropy switch)
   10. Apply radiative cooling
   11. Assign new particle timesteps
   12. Output snapshots / diagnostics
}
```

This corresponds to the Kick-Drift-Kick (KDK) operator-split leapfrog scheme described in Section 10.



## 3. Cosmological Framework

### 3.1 Supported Cosmologies

The code supports two background cosmologies (as well as a non-expanding static box):

**Einstein-de Sitter (EdS):** $\Omega_m = 1$, $\Omega_\Lambda = 0$. The scale factor evolves as

$$a(t) = t^{2/3}, \qquad H(t) = \frac{2}{3t}$$

in internal code time units.

**Flat $\Lambda$CDM:** Default parameters $\Omega_m = 0.3$, $\Omega_\Lambda = 0.7$, $\Omega_b = 0.045$, $h_{100} = 0.7$. The exact solution for the scale factor in a flat universe with matter and cosmological constant is

$$a(t) = \left(\frac{\Omega_m}{\Omega_\Lambda}\right)^{1/3} \sinh^{2/3}\!\left(\frac{3}{2} H_0 \sqrt{\Omega_\Lambda}\; t\right),$$

which the code evaluates analytically for the time-to-scale-factor mapping. The Hubble parameter as a function of scale factor is

$$H(a) = H_0 \sqrt{\Omega_m a^{-3} + \Omega_\Lambda},$$

evaluated every timestep for Hubble-drag and cosmological expansion-cooling terms.

Additional default cosmological parameters: $\sigma_8 = 0.81$, $n_s = 0.96$ (spectral index).

### 3.2 Comoving Coordinates

All positions are stored in comoving coordinates $\mathbf{x} = \mathbf{r}/a$, and velocities are stored as comoving velocities $\mathbf{v} = \dot{\mathbf{x}}$. During each drift step, positions are updated in comoving space. During each kick of duration $\Delta t/2$, the Hubble drag term

$$\left.\dot{\mathbf{v}}\right\vert{}_{\mathrm{drag}} = -2H(a)\,\mathbf{v}$$

is applied multiplicatively to all particle velocities via the discrete factor $(a / a_\mathrm{next})^2$, where $a_\mathrm{next} = a(1 + H\Delta t/2)$.

### 3.3 Internal Unit System and Physical Conversions

Internally, the simulation evolves in dimensionless code units defined by fixing the base domain side length to $1$, the total box mass to $M_\mathrm{tot} = 1$, and the dimensionless Hubble constant to $H_0^\mathrm{code} = \frac{2}{3}\Omega_m^{-1/2}$. Under this normalization, the gravitational constant in code units is fixed to $G = 1/(6\pi)$.

To interface with physical quantities (such as the initial power spectrum, radiative cooling tables, temperature floor, and snapshot metadata), the user specifies the physical box size at $a=1$ via the configuration parameter ($L_\mathrm{box}$). This defines the physical equivalence of $1.0$ code unit:

| Quantity | Physical equivalence of $1.0$ code unit |
|  |  |
| Length | $U_L = L_\mathrm{box}$ (comoving Mpc, i.e., physical Mpc at $a=1$) |
| Time | $U_t = H_0^\mathrm{code} / H_0^\mathrm{phys}$ (Gyr) |
| Mass | $U_M = \rho_\mathrm{crit}\,\Omega_m\,L_\mathrm{box}^3$ ($M_\odot$) |
| Velocity | $U_v = U_L / U_t = 150\, h_{100}\sqrt{\Omega_m}\,(L_\mathrm{box}/\mathrm{Mpc})\;\mathrm{km\,s^{-1}}$ (comoving) |

Physical temperatures are recovered from the comoving code specific internal energy $u_\mathrm{code}$ via

$$T = a^2\, u_\mathrm{code}\, U_v^2\, (\gamma - 1)\,\frac{\mu\, m_p}{k_B},$$

where $\mu = 1.22$ is the mean molecular weight for a neutral primordial H/He plasma (76% H, 24% He), $m_p$ is the proton mass, and $k_B$ is the Boltzmann constant.



## 4. Initial Conditions

### 4.1 Power Spectrum

Initial density perturbations are drawn from a Gaussian random field with power spectrum

$$P(k) = A_s\, k^{n_s}\, T^2(k),$$

where $T(k)$ is the BBKS transfer function (Bardeen et al. 1986)

$$T(q) = \frac{\ln(1+2.34q)}{2.34q}\left[1 + 3.89q + (16.1q)^2 + (5.46q)^3 + (6.71q)^4\right]^{-1/4},$$

with $q = k / (\Gamma\, h\,\mathrm{Mpc}^{-1})$ and shape parameter $\Gamma = \Omega_m h_{100}$. The amplitude $A_s$ is determined by normalizing to $\sigma_8$:

$$\sigma_8^2 = \frac{1}{2\pi^2}\int_0^\infty P(k)\, W^2(kR_8)\, k^2\, dk, \qquad R_8 = 8\,h^{-1}\,\mathrm{Mpc},$$

where $W(x) = 3(\sin x - x\cos x)/x^3$ is the spherical top-hat window function.

### 4.2 Linear Growth Factor

The present-day linear growth suppression factor $g_1(\Omega_{m,0})$ is computed using the Carroll, Press & Turner (1992) fitting formula:

$$g_1(\Omega_m) \approx \frac{5\Omega_m}{2}\left[\Omega_m^{4/7} - \Omega_\Lambda + \left(1+\frac{\Omega_m}{2}\right)\left(1+\frac{\Omega_\Lambda}{70}\right)\right]^{-1}.$$

Because cosmological runs start in the matter-dominated era where $g_1(\Omega_m(a_\mathrm{start})) \approx 1$, the code scales the present-day displacement field back to the starting epoch via $D_1(a_\mathrm{start}) = a_\mathrm{start} / g_1(\Omega_{m,0})$. The logarithmic growth rate used for initial velocities is approximated by the Peebles (1980) formula

$$f \equiv \frac{d\ln D_1}{d\ln a} \approx \Omega_m(a)^{0.55}, \qquad \Omega_m(a) = \frac{\Omega_m}{\Omega_m + \Omega_\Lambda a^3}.$$

### 4.3 Zel'dovich Approximation (1LPT)

Particle positions and velocities are initialized using the Zel'dovich approximation (Zel'dovich 1970), i.e., first-order Lagrangian Perturbation Theory (1LPT). Given a Lagrangian coordinate $\mathbf{q}$, the perturbed comoving position at scale factor $a$ is

$$\mathbf{x}(\mathbf{q}, a) = \mathbf{q} + \mathbf{d}(\mathbf{q}, a), \qquad \mathbf{d}(\mathbf{q}, a) = -D_1(a)\,\nabla_q \Phi(\mathbf{q}),$$

where $\Phi(\mathbf{q})$ is the displacement potential satisfying $\nabla^2\Phi = \delta$. In Fourier space the displacement is

$$\tilde{\mathbf{d}}(\mathbf{k}) = -i\frac{\mathbf{k}}{k^2}\,\tilde{\delta}(\mathbf{k}),$$

and the comoving velocity field is

$$\mathbf{v}(\mathbf{q}, a) = H(a)\, f(a)\, \mathbf{d}(\mathbf{q}, a).$$

A real-space Gaussian white-noise field is generated on the mesh, transformed to Fourier space via a 3D r2c FFT, scaled by $\sqrt{P(k)}$, deconvolved by the CIC interpolation window ($\mathrm{sinc}^{-2}$ per dimension), transformed back to real space, and interpolated to particle positions using CIC. The code optionally supports:

- **Fixed-amplitude ICs**: Fourier amplitudes fixed to the theoretical expectation $\sqrt{P(k)}$, with random phases preserved (Angulo & Pontzen 2016).
- **Phase inversion**: $\tilde{\delta}(\mathbf{k}) \to -\tilde{\delta}(\mathbf{k})$, enabling paired-and-fixed realizations.
- **Glass files**: pre-generated glass particle cubes read in HDF5 format and tiled across the domain to suppress lattice artifacts.



## 5. Gravity Solver

### 5.1 PM+PP Force Split

The total gravitational force on particle $i$ is split into a long-range PM component and a short-range PP component:

$$\mathbf{F}_i = \mathbf{F}_i^\mathrm{PM,long} + \mathbf{F}_i^\mathrm{PP,short}.$$

The split is defined by a scale radius $r_s = N_\mathrm{smooth}\,\Delta x$, where $\Delta x = 1/N_\mathrm{grid}$ is the PM cell size in code units and $N_\mathrm{smooth} = 1.25$ by default. The PP kernel corrects for the under-resolved near-field that the PM mesh misses, acting within a cutoff radius $r_\mathrm{cut} = 4.5\,r_s$.

### 5.2 Particle-Mesh Solver

**Mass assignment.** Cloud-In-Cell (CIC) interpolation assigns each particle's mass to the nearest $2^3 = 8$ cell-centered grid nodes, producing a density field $\rho(\mathbf{x})$ on the PM grid.

**Poisson equation.** The gravitational potential is solved in Fourier space via

$$\hat{\phi}(\mathbf{k}) = -\frac{4\pi G\,\hat{\rho}(\mathbf{k})}{k^2}\, S(\mathbf{k}),$$

where the Green's function includes a smoothing and deconvolution window

$$S(\mathbf{k}) = \frac{\exp(-k^2 r_s^2)}{\mathrm{sinc}^4(k_x\Delta x/2)\,\mathrm{sinc}^4(k_y\Delta x/2)\,\mathrm{sinc}^4(k_z\Delta x/2)},$$

which simultaneously deconvolves the CIC mass assignment and force interpolation ($\mathrm{sinc}^{-4}$) and exponentially damps the long-range force below $r_s$ so the PP correction can take over. The $k = 0$ mode is set to zero (periodic boundaries, zero mean field).

**Force from potential.** The comoving acceleration is obtained by differentiating the potential in real space using a 4-point central-difference stencil with $O(\Delta x^4)$ accuracy, scaled by $a^{-3}$:

$$\mathbf{g}(\mathbf{x}_i) = -\frac{1}{a^3\,\Delta x}\left[\frac{2}{3}\left(\phi_{+1}-\phi_{-1}\right) - \frac{1}{12}\left(\phi_{+2}-\phi_{-2}\right)\right]\hat{\mathbf{e}},$$

applied independently along each axis. The acceleration is then interpolated back to particle positions using the cached CIC weights.

### 5.3 Particle-Particle Short-Range Correction

The PP correction provides the near-field gravity missing from the PM mesh. For each particle pair $(i, j)$ with $r_{ij} < r_\mathrm{cut}$, the short-range correction force is

$$\mathbf{F}_{ij}^\mathrm{PP} = -G_\mathrm{com}\, m_i m_j\, \frac{\partial \phi_\mathrm{soft}(r, \varepsilon)}{\partial r}\, f_\mathrm{taper}(r)\,\hat{\mathbf{r}}_{ij},$$

where $G_\mathrm{com} = G / a^3$, $\partial \phi_\mathrm{soft}/\partial r$ is the cubic-spline softened force derivative (Section 5.4), and the complementary error-function force taper is

$$f_\mathrm{taper}(r) = \mathrm{erfc}\!\left(\frac{r}{2r_s}\right) + \frac{r}{\sqrt{\pi}\,r_s}\exp\!\left(-\frac{r^2}{4r_s^2}\right).$$

This splitting is the exact real-space counterpart to the Gaussian-smoothed PM Green's function, ensuring that the long-range and short-range components sum smoothly to the full softened $1/r^2$ force.

### 5.4 Gravitational Softening

The force is softened to suppress two-body scattering at scales below the mean inter-particle spacing. The base comoving softening length is

$$\varepsilon_\mathrm{base} = f_\varepsilon\,\bar{d},$$

where $\bar{d} = 1 / N_\mathrm{part,1D}$ is the mean inter-particle spacing in code units and $f_\varepsilon = 0.033$. For scale factors $a > a_\mathrm{cap} = 0.3$, the comoving softening is scaled as $\varepsilon(a) = \varepsilon_\mathrm{base}(a_\mathrm{cap}/a)$ so that the physical softening remains constant, following the convention of Springel (2005).

The softened gravitational force uses the $1h$-support cubic-spline kernel derivative (Price & Monaghan 2007):

$$\frac{\partial \phi(r, \varepsilon)}{\partial r} = \frac{1}{\varepsilon^2} \begin{cases} \frac{32}{3}q - \frac{192}{5}q^3 + 32q^4, & 0 \le q < \frac{1}{2}, \\ -\frac{1}{15q^2} + \frac{64}{3}q - 48q^2 + \frac{192}{5}q^3 - \frac{32}{3}q^4, & \frac{1}{2} \le q < 1, \\ \frac{1}{q^2}, & q \ge 1, \end{cases}$$

where $q = r/\varepsilon$.

**Adaptive softening (gas).** When adaptive softening is enabled, each gas particle uses an adaptive softening length $\varepsilon_i = \mathrm{clamp}(h_i, \varepsilon_{\min}, \varepsilon_{\max})$, where $\varepsilon_{\min} = \varepsilon(a)$ and $\varepsilon_{\max} = r_\mathrm{cut}$ (when PM is active). Pairwise forces are symmetrized via $\frac{1}{2}[\phi'(r,\varepsilon_i) + \phi'(r,\varepsilon_j)]$, and the Price & Monaghan (2007) energy-conservation correction term

$$\mathbf{a}_i^\mathrm{corr} = -\frac{G_\mathrm{com}}{2 m_i} \sum_j \left(m_i \xi_i \frac{\partial W(r_{ij}, \varepsilon_i)}{\partial r} + m_j \xi_j \frac{\partial W(r_{ij}, \varepsilon_j)}{\partial r}\right)\hat{\mathbf{r}}_{ij}$$

is added to gas–gas pairs, where $\xi_i = \frac{m_i h_i}{3\rho_i\Omega_i}\sum_j m_j \frac{\partial\phi(r_{ij}, h_i)}{\partial h_i}\,\mathrm{erfc}(r_{ij}/(2r_s))$.



## 6. LBVH Tree

### 6.1 Construction Algorithm

Short-range PP gravity and hydrodynamic neighbor searches use a Linear Bounding Volume Hierarchy (LBVH) built with the Karras (2012) algorithm. Construction proceeds in five stages:

1. **Morton encoding.** Each particle's 3D periodic position is mapped to a 63-bit Morton (Z-order curve) key by interleaving 21-bit fixed-point representations of the three coordinates:

$$z_i = \mathrm{morton3D}(x_i^{(1)}, x_i^{(2)}, x_i^{(3)}) \in [0, 2^{63}).$$

2. **Sorting.** Particles are sorted by Morton key, bringing spatially proximate particles into contiguous memory order. On the CPU this uses a sort function on an index map followed by parallel SoA permutation; on the GPU a parallel bitonic sort on power-of-two padded buffers is used.

3. **Topology construction.** The binary tree topology is constructed in $O(N)$ using the parallel algorithm of Karras (2012): each internal node $k$ determines its range $[i, j]$ of Morton keys via longest-common-prefix queries and assigns left/right child and parent pointers without locks.

4. **Bounding-box aggregation.** Axis-aligned bounding boxes (AABBs) are propagated bottom-up from leaves (padded by $h_i$ for gas trees) to the root. On the CPU this uses an atomic arrival counter per internal node; on the GPU a level-synchronous bottom-up sweep is used.

5. **Center-of-mass aggregation.** Total mass and center-of-mass $\mathbf{r}_\mathrm{cm}$ are accumulated bottom-up alongside the AABBs to support future multipole extensions.

### 6.2 Tree Traversal

During PP gravity and hydrodynamic neighbor searches, the tree traversal uses a fixed-size per-thread private stack (size 128) and an exact periodic minimum-distance AABB test. Nodes whose bounding boxes lie outside the search radius ($r_\mathrm{cut}$ for PP gravity, $h_i$ for hydro) are pruned, and all intersecting leaf particles are evaluated directly.

### 6.3 GPU Variant

When compiled with GPU support, the PP gravity tree-traversal loop offloads to the GPU, using a per-thread private traversal stack. Unified DM+gas arrays and the LBVH remain resident in device memory, and the full LBVH construction (Morton encoding, bitonic sort, Karras topology, and level-synchronous bottom-up aggregation) executes directly on the GPU.



## 7. Hydrodynamics: Meshless Finite Mass

The hydrodynamics solver implements the Meshless Finite Mass (MFM) method of Hopkins (2015). MFM is a quasi-Lagrangian scheme that evolves the fluid equations in conservation form without a fixed mesh, inheriting SPH's Galilean invariance and adaptive resolution while achieving second-order accuracy in smooth flows and accurately capturing shocks via a Riemann solver.

### 7.1 Kernel and Density Estimation

The 3D cubic-spline kernel with $1h$ compact support is

$$W(r, h) = \frac{8}{\pi h^3} \begin{cases} 1 - 6\left(\frac{r}{h}\right)^2 + 6\left(\frac{r}{h}\right)^3, & 0 \leq \frac{r}{h} < \frac{1}{2}, \\ 2\left(1-\frac{r}{h}\right)^3, & \frac{1}{2} \leq \frac{r}{h} < 1, \\ 0, & \frac{r}{h} \geq 1. \end{cases}$$

The smoothing length $h_i$ for particle $i$ is determined by requiring the kernel-weighted neighbor count to match a target $N_\mathrm{ngb} = 32$:

$$N_\mathrm{eff}(h_i) = \frac{4\pi}{3} h_i^3 \sum_j W(r_{ij}, h_i) = N_\mathrm{ngb}.$$

This implicit equation is solved by Newton-Raphson iteration with geometric bisection fallback and a per-thread neighbor cache to minimize tree traversals. The particle density is then $\rho_i = m_i \sum_j W(r_{ij}, h_i)$ (clamped above the density floor $= 10^{-12}$). After updating $h_i$, the gas LBVH bounding boxes are refreshed.

### 7.2 Gradient Estimation

For each particle $i$, the volume-weighted least-squares gradient matrix is constructed as

$$E_{\alpha\beta}^{(i)} = \sum_{j \neq i} V_j \left(x_j^\alpha - x_i^\alpha\right)\left(x_j^\beta - x_i^\beta\right) W(r_{ij}, h_i),$$

where $V_j = m_j/\rho_j$ is the effective particle volume. The inverse $\mathbf{B}_i = \mathbf{E}_i^{-1}$ is computed directly for each $3\times 3$ matrix, and its condition number is evaluated as $\kappa_i = \frac{1}{3}\|\mathbf{E}_i\|_F\|\mathbf{B}_i\|_F$ (an optional Tikhonov-regularized eigendecomposition pseudoinverse is also available). If $\mathbf{E}_i$ is singular or $\kappa_i$ exceeds the maximum allowed condition number (default $5.0$), particle $i$ is flagged as ill-conditioned and its gradients are set to zero. For well-conditioned particles, the gradient of any primitive scalar $\phi \in \{\rho, P, v_x, v_y, v_z\}$ at particle $i$ is

$$(\nabla\phi)_i = \mathbf{B}_i \sum_{j \neq i} V_j (\phi_j - \phi_i)\, \mathbf{x}_{ij}\, W(r_{ij}, h_i),$$

where $\mathbf{x}_{ij} = \mathbf{x}_j - \mathbf{x}_i$.

### 7.3 Slope Limiting

Gradients are limited using a Barth-Jespersen (1989) slope limiter adapted for particle methods. For density and pressure, the limiter coefficient $\alpha_i \in [0, 1]$ bounds the maximum and minimum extrapolated face differences against the neighborhood extrema scaled by $\beta \in [\beta_{\min}, \beta_{\max}] = [1.0, 2.0]$:

$$\beta = \max\!\left(\beta_{\min},\; \beta_{\max}\,\min\!\left(1, \frac{\kappa_\mathrm{crit}}{\kappa_i}\right)\right),$$

providing stronger limiting in poorly-conditioned regions and less dissipation where $\mathbf{E}_i$ is well-conditioned. Velocity gradients are limited isotropically using the magnitude of the extrapolated velocity difference vector.

At each face, a secondary pairwise limiter is applied to $\rho$ and $P$ to guarantee monotonicity across the interface, while velocity reconstruction uses a Galilean- and rotationally-invariant limiter based on projecting velocities onto the face normal $\hat{\mathbf{n}}$.

### 7.4 Face Reconstruction and Area Vector

The effective face between particles $i$ and $j$ is located at the smoothing-length-weighted midpoint

$$\mathbf{x}_f = \mathbf{x}_i + \frac{h_i}{h_i + h_j}\,\mathbf{x}_{ij}.$$

The effective MFM face area vector (pointing from $i$ to $j$) is

$$\mathbf{A}_{ij} = V_i V_j W(r_{ij}, h_i)\,\mathbf{B}_i\mathbf{x}_{ij} + V_i V_j W(r_{ij}, h_j)\,\mathbf{B}_j\mathbf{x}_{ij}.$$

If $\mathbf{A}_{ij}\cdot\mathbf{x}_{ij} < 0$ or either particle is ill-conditioned, the solver falls back to zeroth-order (piecewise-constant) reconstruction using the SPH-like area vector

$$\mathbf{A}_{ij}^\mathrm{SPH} = -\left(V_i^2\, \frac{\partial W_i}{\partial r} + V_j^2\, \frac{\partial W_j}{\partial r}\right) \frac{\mathbf{x}_{ij}}{r_{ij}}.$$

Otherwise, left and right states at the face are obtained by linear extrapolation from $i$ and $j$ to $\mathbf{x}_f$ using the limited gradients and pairwise limiter, clamped above density and pressure floors.

### 7.5 Riemann Solver

#### HLLC Solver (default)

The HLLC solver (Toro et al. 1994) is used at each face in the frame moving at velocity $\mathbf{v}_f = \mathbf{v}_i + \frac{h_i}{h_i+h_j}(\mathbf{v}_j - \mathbf{v}_i)$. The solver first computes Roe-averaged Einfeldt (1988) wave speed bounds:

$$S_L = \min(v_{n,L} - c_L,\; v_\mathrm{Roe} - c_\mathrm{Roe}), \qquad S_R = \max(v_{n,R} + c_R,\; v_\mathrm{Roe} + c_\mathrm{Roe}).$$

If $S_L \ge S_R$ or the resulting star-region pressure is non-positive ($P^* \le 0$), the solver falls back to simple Davis bounds $S_L = \min(v_{n,L} - c_L, v_{n,R} - c_R)$ and $S_R = \max(v_{n,L} + c_L, v_{n,R} + c_R)$. In the subsonic regime ($S_L < 0 < S_R$), the contact wave speed and interface pressure are

$$S^* = \frac{P_R - P_L + \rho_L v_{n,L}(S_L - v_{n,L}) - \rho_R v_{n,R}(S_R - v_{n,R})}{\rho_L(S_L - v_{n,L}) - \rho_R(S_R - v_{n,R})},$$

$$P^* = P_L + \rho_L(S_L - v_{n,L})(S^* - v_{n,L}).$$

#### Exact Riemann Solver (optional)

An exact Riemann solver following Toro (2009, Chapter 4) is available as a runtime configuration option. It iterates Newton-Raphson on the pressure equation

$$f(P^*) = f_L(P^*, \mathbf{U}_L) + f_R(P^*, \mathbf{U}_R) + (v_{n,R} - v_{n,L}) = 0,$$

handling shock, rarefaction, and vacuum branches separately to a relative tolerance of $10^{-6}$ (maximum 100 iterations). If either solver fails on reconstructed states, the code retries the exact solver using first-order (cell-center) states; if that also fails, a zero-flux fallback is applied and logged.

### 7.6 Conservation Equations

Because MFM has zero mass flux across the moving face, each face interaction yields a momentum flux $\mathbf{F}_{ij}^\mathrm{mom} = P^*\mathbf{A}_{ij}$ and a total energy flux $\dot{E}_{ij} = P^*\,\mathbf{v}^*_\mathrm{lab}\cdot\mathbf{A}_{ij}$, where $\mathbf{v}^*_\mathrm{lab} = \mathbf{v}_f + S^*\hat{\mathbf{n}}$ is the lab-frame contact velocity. The resulting rates of change for particle $i$ are

$$\frac{d(m_i \mathbf{v}_i)}{dt} = -\sum_j P^*_{ij}\,\mathbf{A}_{ij} + m_i \mathbf{a}_i^\mathrm{grav},$$

$$\frac{d(m_i E_i)}{dt} = -\sum_j P^*_{ij}\,\mathbf{v}^*_{\mathrm{lab},ij}\cdot\mathbf{A}_{ij} + m_i \mathbf{a}_i^\mathrm{grav}\cdot\mathbf{v}_i,$$

with the corresponding internal energy work rate $\dot{u}_i = -\frac{1}{m_i}\sum_j P^*_{ij}(\mathbf{v}^*_{\mathrm{lab},ij} - \mathbf{v}_i)\cdot\mathbf{A}_{ij}$. By construction, $\mathbf{A}_{ij} = -\mathbf{A}_{ji}$, guaranteeing exact antisymmetry. Pairwise impulses are accumulated using the minimum interacting timestep $\Delta t_\mathrm{int} = \min(\Delta t_i, \Delta t_j)$ to maintain exact momentum and energy conservation under individual particle timesteps.

### 7.7 Equation of State

The gas obeys an ideal-gas equation of state:

$$P_i = (\gamma - 1)\,\rho_i\, u_i, \qquad \gamma = 5/3.$$

The sound speed is $c_i = \sqrt{\gamma P_i/\rho_i}$. The entropic function $A_i = P_i / \rho_i^\gamma$ is stored alongside $u_i$ and $E_i$ as an auxiliary thermodynamic variable.



## 8. Thermodynamics and Dual-Energy Formalism

### 8.1 Energy Arbitration

In cosmological simulations spanning many orders of magnitude in temperature, the internal energy $u$ can be a tiny fraction of the total energy $E$, leading to catastrophic cancellation when computing $u = E - \frac{1}{2}v^2$ in cold, high-Mach-number flows. Following Bryan et al. (1995) and Steinborn et al. (2016), we implement a Dual Energy Formalism (DEF) with an entropy switch that arbitrates each kick step using the explicitly integrated internal energy $u_\mathrm{raw}$, the specific kinetic energy $KE_\mathrm{spec}$, and the local gravitational work scale $\Delta E_\mathrm{grav} = |\mathbf{a}_\mathrm{grav}| h_i$:

1. **Conservative regime (hot/shocked):** If $u_\mathrm{raw} > \alpha_\mathrm{DEF}\,(KE_\mathrm{spec} + \Delta E_\mathrm{grav})$ (default $10^{-3}$), set $u_i = \max(E_i - KE_\mathrm{spec}, 0)$.
2. **Entropy fallback (cold/smooth):** Otherwise, if $u_\mathrm{raw} < \alpha_\mathrm{kin}\,(\max_{j\in\mathcal{N}(i)} KE_{\mathrm{rel},ij} + u_\mathrm{raw})$ or $u_\mathrm{raw} < \alpha_\mathrm{grav}\,\Delta E_\mathrm{grav}$ (both default $10^{-3}$), recover $u_i$ adiabatically from the stored entropic function $A_i$:

$$u_i = \frac{A_i\,\rho_i^{\gamma - 1}}{\gamma - 1},$$

and resynchronize $E_i = u_i + KE_\mathrm{spec}$.

3. **DEF regime (intermediate):** Otherwise, retain the explicitly integrated $u_\mathrm{raw}$ and resynchronize $E_i = u_i + KE_\mathrm{spec}$.

### 8.2 Cosmological Adiabatic Cooling

Gas in an expanding universe undergoes adiabatic $P\,dV$ cooling. Because the code evolves comoving specific internal energy $u_\mathrm{code} = u_\mathrm{phys}/a^2$, each kick applies the exact discrete expansion factor

$$u_i \leftarrow u_i \left(\frac{a}{a_\mathrm{next}}\right)^{3\gamma - 1},$$

which combines physical adiabatic cooling ($a^{-3(\gamma-1)}$) and comoving energy scaling ($a^{-2}$).

### 8.3 Temperature Floor

A constant physical temperature floor is enforced at $T_\mathrm{floor} = 2.73\;\mathrm{K}$ (configurable), converted to comoving code units at each scale factor $a$ via $u_\mathrm{floor} = T_\mathrm{floor} / (a^2 f_{u\to T})$.



## 9. Radiative Cooling

### 9.1 Tabulated UV-Background Cooling

The primary cooling model reads pre-computed Grackle-format HDF5 tables (Smith et al. 2017) for the Haardt & Madau (2012; HM2012) UV background, storing primordial and metal cooling and heating rates on a 3D grid over hydrogen number density $\log_{10} n_H\;(\mathrm{cm}^{-3})$, redshift $z$, and temperature $\log_{10} T\;(\mathrm{K})$.

At runtime, the rates are evaluated by trilinear interpolation in $(z,\,\log_{10} n_H,\,\log_{10} T)$ and combined according to local metallicity as

$$\Lambda_\mathrm{tot} = \Lambda_\mathrm{prim} + \frac{Z}{Z_\odot}\,\Lambda_\mathrm{metal}, \qquad \Gamma_\mathrm{tot} = \Gamma_\mathrm{prim} + \frac{Z}{Z_\odot}\,\Gamma_\mathrm{metal},$$

with $Z_\odot = 0.02$.

### 9.2 Implicit Cooling Integration

In physical CGS units, the specific internal energy evolves due to radiative cooling and photoheating as

$$\frac{du}{dt} = \frac{n_H^2}{\rho}\left[\Gamma_\mathrm{tot}(z, n_H, T) - \Lambda_\mathrm{tot}(z, n_H, T)\right], \qquad n_H = \frac{X\rho}{m_p}\;(X = 0.76).$$

In comoving code units, this equation is integrated implicitly over each active particle's timestep $\Delta t$ using a bisection solver with relative tolerance $10^{-5}$ and a maximum of 50 iterations. This guarantees unconditional stability across the steep cooling curve near $T \sim 10^4\;\mathrm{K}$.

### 9.3 Analytic Fallback

When tabulated cooling data are unavailable, an analytic Bremsstrahlung cooling function is used as a fallback:

$$\Lambda_\mathrm{ff} = 1.4 \times 10^{-27}\, n^2\, \sqrt{T}\quad [\mathrm{erg\,s^{-1}\,cm^{-3}}], \qquad n = \frac{\rho}{\mu m_p},$$

operating above the cooling cutoff ($10^4\;\mathrm{K}$). Cooling is applied after Kick 2.



## 10. Time Integration

### 10.1 Kick-Drift-Kick Leapfrog

The code uses a second-order Kick-Drift-Kick (KDK) leapfrog integrator. For a step of size $\Delta t$:

1. **Kick 1:** Advances DM and gas velocities by $\Delta t/2$ using gravitational accelerations and applies Hubble drag and expansion cooling (saving pre-kick gas snapshots for potential wake-up rollbacks).
2. **Drift:** Advances all particle comoving positions by $\Delta t$ with periodic wrapping.
3. **Force & Flux Evaluation:** Updates cosmology, rebuilds LBVH trees, computes PM gravity, solves gas densities and smoothing lengths, predicts gas velocities and adiabatic pressures to the current global time, computes MFM gradients and Riemann impulses, and evaluates PP gravity.
4. **Kick 2:** Advances DM and gas velocities by $\Delta t/2$ using the new gravitational accelerations, applies the accumulated hydrodynamic impulses to gas particles, applies Hubble drag and expansion cooling, and performs DEF/entropy arbitration.
5. **Cooling:** Applies implicit radiative cooling to active gas particles.

### 10.2 Individual Particle Timesteps

An optional hierarchical block timestep scheme (Springel 2005) assigns each particle a power-of-two integer timestep from a global hierarchy with $2^{30}$ ticks spanning the base dynamical step.

The timestep limiter of Saitoh & Makino (2009) is enforced: if neighboring gas particles differ in timestep by more than a factor of 4, inactive neighbors are flagged for forced wake-up. On the next sync, woken particles roll back to their saved pre-sleep state and perform a catch-up kick over the exact elapsed sleep interval.

### 10.3 Timestep Criteria

Each particle's ideal timestep is the minimum of the applicable physical criteria (capped by the maximum timestep, snapped to the nearest power of 2 below the maximum dynamical timestep factor $\times G^{-1/2}$):

**Gravitational:**
$$\Delta t_\mathrm{grav} = \sqrt{\frac{2\,\eta\,\varepsilon_i}{|\mathbf{a}_\mathrm{grav}|}}, \qquad \eta = 0.01,$$

where $\varepsilon_i$ is the particle's softening length.

**CFL (hydro):**
$$\Delta t_\mathrm{CFL} = C_\mathrm{CFL}\,\frac{h_i}{v_{\mathrm{sig},\max,i}}, \qquad C_\mathrm{CFL} = 0.3,$$

where the pairwise signal velocity is

$$v_{\mathrm{sig},ij} = c_i + c_j - \min\!\left(0,\; \hat{\mathbf{r}}_{ij}\cdot(\mathbf{v}_j - \mathbf{v}_i)\right),$$

and $v_{\mathrm{sig},\max,i}$ is the maximum over all interacting neighbors.

**Expansion:**
$$\Delta t_\mathrm{exp} \leq \frac{0.01}{H(a)},$$

preventing the timestep from exceeding $1\%$ of the Hubble time. Because radiative cooling is integrated implicitly, no explicit cooling timestep limiter is required.



## 11. Code Tests

### 11.1 Sod Shock Tube

The standard Sod (1978) shock-tube test consists of two uniform gas states separated by a discontinuity at $x = L_x/2$. Left state: $(\rho_L, v_L, P_L) = (1, 0, 1)$; right state: $(\rho_R, v_R, P_R) = (0.125, 0, 0.1)$. Particles are initialized using tiled glass cubes (with static boundary zones at $x < 0.05 L_x$ and $x > 0.95 L_x$) and evolved to $t = 1.8$.

### 11.2 Sedov–Taylor Blast Wave

The Sedov–Taylor test (Sedov 1959; Taylor 1950) initializes a point explosion of energy $E_\mathrm{blast} = 1$ in a uniform background ($\rho_0 = 1$, $P_0 = 10^{-6}$). The analytic self-similar blast wave radius is

$$R_s(t) = \xi_0 \left(\frac{E_\mathrm{blast}}{\rho_0}\right)^{1/5} t^{2/5},$$

where $\xi_0 \approx 1.033$ for $\gamma = 5/3$. We inject $E_\mathrm{blast}$ as thermal energy into the central particles (Gaussian-weighted over the 15 nearest glass particles) and verify the shock position, peak density ($\rho_\mathrm{max}/\rho_0 = 4$ for $\gamma = 5/3$), and velocity profile against the analytic solution.

### 11.3 Sound Wave Propagation

The linear sound wave test (Groth et al. 2023, OG3 test) propagates a plane wave with amplitude $\delta\rho/\rho_0 = 10^{-4}$ and sound speed $c_s = 2/3$ in a periodic box. The measured phase speed and damping rate are compared with the analytic solution, confirming second-order convergence ($L_1$ error $\propto h^2$).

### 11.4 Kelvin–Helmholtz Instability

We reproduce the McNally, Lyra & Passy (2012) Kelvin–Helmholtz test in a periodic 3D box with smoothed density ($\rho_1 = 1$, $\rho_2 = 2$) and shear velocity ($v_1 = -0.5$, $v_2 = 0.5$) transitions of width $\delta_y = 0.025$ and sinusoidal vertical velocity perturbation $v_y = 0.01 \sin(4\pi x/L_x)$. The linear growth rate is compared against the analytic expression

$$\omega_\mathrm{KH} = k\,|\Delta v|\,\frac{\sqrt{\rho_1\rho_2}}{\rho_1 + \rho_2},$$

and the nonlinear vortex roll-up morphology is compared with the reference solution of Hopkins (2015).

### 11.5 Evrard Collapse

The Evrard (1988) collapse test initializes a cold sphere of gas (radius $R = 0.35$, mass $M = \Omega_b/\Omega_m$ in code units) with a $\rho(r) \propto 1/r$ radial density profile (generated by mapping a uniform sphere via $r_\mathrm{target} = R\,(r_u/R)^{3/2}$) and specific internal energy $u = 0.05\,GM/R$. The sphere collapses, bounces, and settles into hydrostatic equilibrium. Radial velocity and entropy profiles are compared against the standard reference solution of Hernquist & Katz (1989).

### 11.6 Zel'dovich Pancake

The 1D cosmological Zel'dovich (1970) pancake test evolves a sinusoidal perturbation collapsing at caustic redshift $z_c = 1$ in a periodic Einstein-de Sitter box from initial scale factor $a_i$. The analytic shell-crossing epoch is reproduced within the limits of the 1D resolution, the post-collapse density profile follows the expected caustic cusp, and the pre-collapse velocity field agrees with the linear theory prediction $v(q, a) = -a H(a) f(a) D_1(a)\,\partial\Phi/\partial q$.



## 12. Performance

### 12.1 OpenMP Scaling

All major kernels — PP gravity tree-traversal, density iteration, gradient computation, Riemann flux loops, and leapfrog updates — are parallelized with OpenMP. Dynamic scheduling balances load across irregular neighbor counts in clustered regions, while static scheduling with SIMD vectorization minimizes overhead in regular loops.

### 12.2 GPU Acceleration

When compiled with GPU support, the PP gravity tree traversal offloads to the GPU, avoiding vendor-specific extensions and enabling portability across different GPU backends. Unified particle arrays and the LBVH remain resident on the device, executing the following GPU kernels:

- Bitonic sort for Morton keys: $O(N\log^2 N)$, fully parallel on padded power-of-two buffers.
- Karras topology build: $O(N)$, embarrassingly parallel per internal node.
- Bottom-up bounding-box and mass sweep: level-synchronous over tree depth.
- Tree-traversal PP force: one GPU thread per active particle, using a per-thread private stack.

The GPU PP gravity pipeline runs concurrently with the CPU MFM hydrodynamics step.

### 12.3 LBVH Cache Efficiency

Sorting the SoA particle arrays into Morton order each step ensures spatial locality: nearby particles in 3D reside nearby in memory, improving cache hit rates during both LBVH neighbor searches and MFM pairwise flux evaluations.



## 13. Data Input/Output

### 13.1 HDF5 Snapshots

Snapshots are written in a custom HDF5 schema using Single-Writer Multiple-Reader (SWMR) mode, allowing concurrent inspection during a running simulation. Each snapshot file contains five groups:

| Group | Contents |
| --- | --- |
| `Header` | Scale factor, simulation time, current softening length and other parameters |
| `Config` | All runtime configuration parameters, domain dimensions, cosmological parameters, and internal $G$ |
| `Units` | Conversion factors to physical Mpc, Gyr, $\mathrm{km\,s^{-1}}$, $M_\odot$, CGS, and $u \to T$ |
| `Particles` | DM positions, velocities, accelerations, masses, and cumulative work attributes |
| `Gas` | Gas positions, velocities, accelerations, masses, $h$, $\rho$, $P$, $u$, entropy, total energy, metallicity, temperature, and diagnostic fields |

Snapshots are written every $\Delta a = 0.005$ by default (configurable).

### 13.2 Diagnostic Output

A CSV log file is written at regular intervals (30 seconds by default) including diagnostic metrics.

### 13.3 Configuration

Runtime parameters are read from an INI-format configuration file. Key parameters and their defaults are summarized in Table 1.

**Table 1.** Key runtime configuration parameters and defaults.

| Section | Parameter | Default | Description |
| --- | --- | --- | --- |
| `[domain]` | `box_size_mpc` | `96.0` | Comoving box side length (Mpc) |
| `[domain]` | `mesh_size_1d` | `64` | PM grid dimension $N_\mathrm{grid}$ |
| `[domain]` | `num_dm_particles_1d` | `64` | DM particles per dimension |
| `[domain]` | `num_gas_particles_1d` | `64` | Gas particles per dimension |
| `[cosmology]` | `omega_m` | `0.3` | Matter density parameter |
| `[cosmology]` | `omega_baryon` | `0.045` | Baryon density parameter |
| `[cosmology]` | `omega_lambda` | `0.7` | Dark energy density parameter |
| `[cosmology]` | `hubble_h` | `0.7` | Dimensionless Hubble constant |
| `[cosmology]` | `sigma_8` | `0.81` | Power spectrum normalization |
| `[cosmology]` | `spectral_index` | `0.96` | Spectral index $n_s$ |
| `[time]` | `a_start` | `0.02` | Starting scale factor |
| `[time]` | `a_end` | `1.0` | Final scale factor |
| `[gravity]` | `comoving_softening_factor` | `0.033` | Softening as fraction of $\bar{d}$ |
| `[gravity]` | `pm_smoothing_cells` | `1.25` | PM smoothing in cells |
| `[gravity]` | `cutoff_radius_factor` | `4.5` | PP cutoff in units of $r_s$ |
| `[hydro]` | `target_neighbors` | `32.0` | Target neighbor count $N_\mathrm{ngb}$ |
| `[time]` | `hydro_courant_factor` | `0.3` | CFL coefficient |
| `[time]` | `gravity_accuracy_eta` | `0.01` | Gravity timestep accuracy $\eta$ |
| `[thermodynamics]` | `gamma` | `5/3` | Adiabatic index |
| `[thermodynamics]` | `dual_energy_threshold` | `0.001` | DEF truncation threshold $\alpha_\mathrm{DEF}$ |
| `[output]` | `save_hdf5_every_delta_a` | `0.005` | Snapshot cadence |



## 14. Summary & Future Work

We have presented a new cosmological hydrodynamics code combining a PM+PP gravity solver accelerated by an LBVH tree with the MFM hydrodynamics scheme. Future work may include multipole approximations in the LBVH traversal, second-order LPT (2LPT) initial conditions, star formation and stellar feedback subgrid models and MHD extensions.



## References

- Angulo, R. E. & Pontzen, A. 2016, MNRAS Letters, 462, L1
- Bardeen, J. M., Bond, J. R., Kaiser, N. & Szalay, A. S. 1986, ApJ, 304, 15 (BBKS)
- Barth, T. J. & Jespersen, D. C. 1989, AIAA Paper 89-0366
- Bryan, G. L. et al. 1995, Computer Physics Communications, 89, 149
- Bryan, G. L. et al. 2014, ApJS, 211, 19 (Enzo)
- Carroll, S. M., Press, W. H. & Turner, E. L. 1992, ARA&A, 30, 499
- Einfeldt, B. 1988, SIAM J. Numer. Anal., 25, 294
- Evrard, A. E. 1988, MNRAS, 235, 911
- Frigo, M. & Johnson, S. G. 2005, Proc. IEEE, 93, 216 (FFTW)
- Groth, F. et al. 2023, A&A, 678, A10 (OG3 test)
- Guennebaud, G. et al. 2010, Eigen v3, http://eigen.tuxfamily.org
- Haardt, F. & Madau, P. 2012, ApJ, 746, 125 (HM2012)
- Hernquist, L. & Katz, N. 1989, ApJS, 70, 419 (TreeSPH)
- Hopkins, P. F. 2015, MNRAS, 450, 53 (GIZMO/MFM)
- Karras, T. 2012, Proc. High Performance Graphics (HPG 2012), 29–37
- Lauterbach, C. et al. 2009, Computer Graphics Forum, 28, 1427
- McNally, C. P., Lyra, W. & Passy, J.-C. 2012, ApJS, 201, 18
- Monaghan, J. J. 1992, ARA&A, 30, 543
- Peebles, P. J. E. 1980, The Large-Scale Structure of the Universe, Princeton Univ. Press
- Potter, D., Stadel, J. & Teyssier, R. 2017, Comput. Astrophys. Cosmol., 4, 2 (pkdgrav3)
- Price, D. J. 2012, J. Comput. Phys., 231, 759
- Price, D. J. & Monaghan, J. J. 2007, MNRAS, 374, 1347
- Saitoh, T. R. & Makino, J. 2009, ApJL, 697, L99
- Schaller, M. et al. 2024, MNRAS, 530, 2378 (SWIFT)
- Sedov, L. I. 1959, Similarity and Dimensional Methods in Mechanics, Academic Press
- Smith, B. D. et al. 2017, MNRAS, 466, 2217 (Grackle)
- Sod, G. A. 1978, J. Comput. Phys., 27, 1
- Springel, V. 2005, MNRAS, 364, 1105 (Gadget-2)
- Springel, V. 2010, MNRAS, 401, 791 (AREPO)
- Springel, V. et al. 2021, MNRAS, 506, 2871 (Gadget-4)
- Steinborn, L. K. et al. 2016, MNRAS, 458, 1540
- Taylor, G. I. 1950, Proc. Roy. Soc. London A, 201, 175
- Teyssier, R. 2002, A&A, 385, 337 (RAMSES)
- Toro, E. F. 2009, Riemann Solvers and Numerical Methods for Fluid Dynamics, 3rd ed., Springer
- Toro, E. F., Spruce, M. & Speares, W. 1994, Shock Waves, 4, 25 (HLLC)
- Zel'dovich, Ya. B. 1970, A&A, 5, 84
