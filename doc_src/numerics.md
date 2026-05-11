# Numerical Method

The LowMachReact-Hex solver uses a collocated finite-volume method on structured hexahedral meshes to solve the incompressible Navier-Stokes and species transport equations.

## Governing Equations

The solver targets the constant-density incompressible formulation:

$$ \nabla \cdot \mathbf{u} = 0 $$

$$ \frac{\partial \mathbf{u}}{\partial t} + \nabla \cdot (\mathbf{u} \mathbf{u}) = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{f} $$

Transport of $N$ passive species is solved as:

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

## Species Transport

*   **Stability:** Species are advanced using explicit upwind advection to ensure boundedness ($0 \le Y_k \le 1$).
*   **Conservation:** After the explicit update, mass fractions are renormalized such that $\sum Y_k = 1$ to maintain consistency.
*   **Diffusivity:** Can be specified as a constant per-species or evaluated dynamically through the **Cantera** bridge.

## Conservation Diagnostics

The solver monitors several diagnostics to ensure stability and convergence:
*   **Max/RMS Divergence:** Monitors the residual of the continuity equation.
*   **Net Boundary Flux:** Verifies global mass conservation.
*   **Kinetic Energy:** Tracks the physical energy evolution of the flow.
