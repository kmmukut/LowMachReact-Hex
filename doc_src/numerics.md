# Numerical Method

The LowMachReact-Hex solver uses a collocated finite-volume method on structured hexahedral meshes to solve the incompressible Navier-Stokes and species transport equations.

## Governing Equations

The solver targets the constant-density incompressible formulation:

$$ \nabla \cdot \mathbf{u} = 0 $$

$$ \frac{\partial \mathbf{u}}{\partial t} + \nabla \cdot (\mathbf{u} \mathbf{u}) = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{f} $$

Transport of $N$ passive species is solved as:
(using a **Scale-on-Demand** architecture supporting $N \ge 0$)

$$ \frac{\partial Y_k}{\partial t} + \nabla \cdot (\mathbf{u} Y_k) = \nabla \cdot (D_k \nabla Y_k) \quad k=1...N $$

## Time Integration

A fractional-step projection method is used for pressure-velocity coupling.

1.  **Momentum Predictor:**
    An intermediate velocity $\mathbf{u}^*$ is computed using an explicit Adams-Bashforth 2 (AB2) scheme:
    $$ \frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t} = \frac{3}{2} \mathbf{R}^n - \frac{1}{2} \mathbf{R}^{n-1} $$
    where $\mathbf{R}$ represents the sum of convective, diffusive, and body force terms.

2.  **Pressure Poisson Solve:**
    The pressure correction $\phi$ is found by solving:
    $$ \nabla^2 \phi = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^* $$
    This equation is solved using a matrix-free **Conjugate Gradient (CG)** method.

3.  **Correction Step:**
    Finally, the velocity and pressure are updated:
    $$ \mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla \phi $$
    $$ p^{n+1} = p^n + \phi $$

### Time Stepping and CFL Control

The solver natively supports both fixed-step and dynamically scaled time integration.
For unstructured/hexahedral grids, the physical cell CFL number is rigorously computed using the absolute sum of conservative boundary fluxes:
$$ CFL_c = \frac{\Delta t}{V_c} \left( \frac{1}{2} \sum_{f} |F_f| \right) $$
By tracking the maximum outward flux ratio across the domain, the solver dynamically scales $\Delta t$ at the beginning of each step to maintain a target $CFL_{max}$, perfectly balancing stability and computational efficiency.

## Finite Volume Discretization

*   **Collocated Layout:** All variables are stored at cell centers.
*   **Convection:** Supports both **1st-order Upwind** and **2nd-order Central** differencing.
*   **Diffusion:** Second-order central difference.
*   **Face Fluxes:** Conservative corrected face fluxes are used to maintain mass conservation and avoid checkerboard oscillations, following the Rhie-Chow principle.
*   **Interpolation:** Distance-weighted interpolation is used for face values to maintain accuracy on non-uniform meshes:
    $$ \psi_f = w \psi_L + (1-w) \psi_R $$
    where $w$ is the geometric weighting factor based on the distance from cell centers to the face.

## Pressure Solver Details

*   **Matrix-Free:** The Laplacian operator is applied directly without storing the full sparse matrix.
*   **Caching:** To optimize performance, neighbor IDs and geometric coefficients ($A_{face}/d_{normal}$) are cached during initialization.
*   **Null-Space Removal:** For purely Neumann systems (e.g., closed cavity), the pressure is pinned at cell 1 to ensure a unique solution.

## Numerical Stability and Convection Schemes

The choice of convection scheme is critical for the stability of the solver, particularly as the flow accelerates or the viscosity decreases.

### Central vs. Upwind Schemes

*   **Central Difference (2nd Order):** Theoretically more accurate but numerically non-dissipative. It is susceptible to "wiggles" or checkerboard oscillations when the local mesh resolution is insufficient.
*   **Upwind Difference (1st Order):** Highly stable and naturally dissipative, but introduces "numerical diffusion" which can smear out sharp gradients.

### The Grid Peclet Number ($Pe$)

The stability of the `central` scheme is governed by the local **Grid Peclet Number**:
$$ Pe_{grid} = \frac{|\mathbf{u}| \Delta x}{\nu} $$
For the central scheme to remain stable and avoid unphysical oscillations, the cell-centered discretization generally requires $Pe_{grid} \le 2$. 

### Circumventing Instability

When a simulation "blows up" using the central scheme (as observed in the `open_channel` case once velocities increased), several strategies can be employed:

1.  **Upwind Convection:** Switching `convection_scheme = "upwind"` provides the necessary numerical dissipation to stabilize the solver at high Peclet numbers.
2.  **Mesh Refinement:** Reducing $\Delta x$ lowers the local Peclet number, potentially allowing the central scheme to remain stable.
3.  **Artificial Dissipation:** Adding a small amount of numerical viscosity (or using flux-limited schemes like TVD/WENO) can provide stability without the full accuracy loss of 1st-order upwind.
4.  **Time-Step Control:** While CFL control handles the temporal stability, spatial stability (like the Peclet limit) is a property of the discretization itself.

## Species Transport

*   **Stability:** Species are advanced using explicit upwind advection to ensure boundedness ($0 \le Y_k \le 1$).
*   **Conservation:** Diffusive fluxes are corrected using a **Correction Velocity** (or diffusive flux correction) to ensure that the net mass flux summed over all species is identically zero ($\sum J_{diff, k} = 0$). After the explicit update, mass fractions are renormalized such that $\sum Y_k = 1$ to maintain consistency and eliminate minor truncation errors.
*   **Diffusivity:** Can be specified as a constant per-species or evaluated dynamically through the **Cantera** bridge.

## Conservation and Performance Diagnostics

The solver monitors several diagnostics to ensure stability, convergence, and performance:
*   **Max/RMS Divergence:** Monitors the residual of the continuity equation.
*   **Net Boundary Flux:** Verifies global mass conservation.
*   **Kinetic Energy:** Tracks the physical energy evolution of the flow.
*   **CFL Number:** Tracks the active global maximum CFL over the domain.
*   **Performance Profiling:** Hierarchical timing of kernels (e.g., `Pressure_Solve`, `Flow_Transport`, `MPI_Communication`) to identify computational bottlenecks.
