title: Developer Guide

# DEVELOPER_GUIDE.md (LowMachReact-Hex)

## Project summary

This repository is a Fortran 2008 MPI finite-volume solver (LowMachReact-Hex) for laminar incompressible / low-Mach-style flow on hexahedral meshes.

## Current Baseline Architecture

The solver currently simulates low-Mach flow using a fractional-step projection method:

*   **Momentum Predictor**: Explicit AB2/Forward-Euler advection and diffusion.
*   **Pressure Poisson Solve**: Matrix-free Conjugate Gradient using distance-weighted face interpolation.
*   **Scale-on-Demand Species**: Dynamic, memory-safe passive multi-species transport with diffusive-flux correction.
*   **Passive Sensible-Enthalpy Energy**: Optional transport of sensible enthalpy `h`; temperature is recovered from `h`, `Y`, and `p0`.
*   **Cantera Thermo Sync**: Optional Cantera recovery of `T`, `cp`, `lambda`, and diagnostic `rho_thermo` from `h,Y,p0`, using a combined sync path where possible.
*   **Open Boundary Support**: Handles Dirichlet velocity inlets, Neumann pressure outlets, fixed species, and fixed-temperature boundaries.
*   **Numerical Robustness**: Corrected flux interpolation for nonuniform hexahedral grids.

> [!NOTE]
> **Pressure Pinning:** In purely periodic or closed (wall-bounded) domains, the pressure matrix has a null space (constant pressure shift). The solver automatically identifies these cases and pins the pressure at cell 1 to ensure convergence. In cases with at least one Dirichlet pressure boundary, pinning is disabled.

Current baseline:

- Cell-centered finite-volume flow solver.
- Projection method with corrected conservative face fluxes.
- Conservative face-flux divergence diagnostic.
- Scale-on-Demand Species Transport:
  - **`mod_species`**: Manages the transport of passive scalars (\(Y_k\)). Supports **Scale-on-Demand** architecture with dynamic allocation for 0 to 256+ species. Implements a **Correction Velocity** (diffusive flux correction) to ensure strict mass conservation when using different species diffusivities (\(D_k\)).
- **`mod_energy`**: Owns sensible-enthalpy transport, temperature recovery, `qrad`, `cp`, thermal conductivity, and diagnostic `rho_thermo`. The transported energy state is `h`, not `T`.
- **`mod_transport_properties`**: Abstracts transport-property evaluation. It provides a bridge to the **Cantera 3.x C++ API** for dynamic evaluation of viscosity and species diffusivity. Flow/projection density remains the configured constant density in the current solver.
- **Cantera thermo bridge**: Supports sensible enthalpy `h(T,Y,p0)`, temperature recovery `T(h,Y,p0)`, and combined thermo sync `(T, cp, lambda, rho_thermo) = sync(h,Y,p0)`.
- **`mod_profiler`**: A hierarchical performance profiling module used to track execution time for critical kernels, Cantera thermo sync, output, and MPI communication. It provides a terminal summary at the end of each simulation.
- **`mod_bc`**: A unified boundary condition manager that supports field-specific types (Velocity, Pressure, Species) for every patch.
- VTU/PVD output working.
- Field-specific boundary-condition split:
  - legacy `patch_type`
  - `patch_velocity_type`
  - `patch_pressure_type`
  - `patch_species_type`
  - `patch_temperature_type`
- MPI model:
  - full mesh replicated on all ranks
  - flow computation decomposed by owned cell ranges
- Current validation targets:
  - lid-driven cavity (Pure, Passive Species, and Cantera modes)
  - periodic body-force channel flow
- Long-term extensions:
  - variable-density low-Mach coupling
  - species-diffusion enthalpy correction
  - reaction source terms and heat release
  - radiation physics with wavenumber-domain decomposition coupled through `qrad`

The current solver should be described as a **laminar incompressible / low-Mach-style constant-density finite-volume solver**, not as a full reacting low-Mach solver yet.

---

## Energy and Cantera thermo rules

Current energy/thermo rules:

```text
h is the transported thermodynamic state
T is recovered from h, Y, and p0
p0 = background_press
rho_flow = params%rho
rho_thermo = diagnostic only
thermo_update_interval = 1
```

When species transport updates composition, preserve transported enthalpy and recover the new temperature:

```text
T_new = T(h_transported, Y_new, p0)
```

Do not preserve old temperature by recomputing:

```text
h_new = h(T_old, Y_new, p0)
```

The energy path may use a combined Cantera thermo sync:

```text
(T, cp, lambda, rho_thermo) = sync(h, Y, p0)
```

This sync must preserve `h`.

Current missing physics:

```text
variable-density low-Mach divergence constraint
reactions and reaction heat release
species-diffusion enthalpy correction: -div(sum_k h_k J_k)
external radiation physics
```

---

## Branching strategy

- `main`: Stable releases and production-ready code.
- `dev`: Primary development branch. All feature development and bug fixes start here.

---

## Non-negotiable rules

- Keep Fortran 2008 compatibility unless explicitly asked otherwise.
- Preserve MPI correctness for both `mpirun -np 1` and multi-rank runs.
- Do not introduce `-ffast-math`.
- Do not remove backward compatibility for legacy `patch_type`.
- Do not change file formats, case input syntax, or output conventions without updating documentation and example case files.
- Do not directly call MPI inside future radiation kernel code.
- Do not put Cantera calls directly inside the flow solver or species transport loops.
- Keep physics changes small and testable.
- Preserve working lid-driven cavity and periodic body-force channel cases.
- Prefer validation-quality incremental changes over large rewrites.
- If a change affects physics, output, MPI, mesh format, case files, or file formats, document the change.
- Public modules, public types, and public procedures must be documented with FORD-compatible comments.
- If a new module is added, update the Makefile dependency/order.
- If a new feature is added, include at least one smoke-test case or documented validation plan.

---

The current flow solver represents incompressible, constant-density, laminar Navier-Stokes with species transport:

$$
\nabla \cdot \mathbf{u} = 0
$$

$$
\frac{\partial \mathbf{u}}{\partial t} + \nabla \cdot (\mathbf{u} \mathbf{u}) = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{f}_{body}
$$

$$
\frac{\partial Y_k}{\partial t} + \nabla \cdot (\mathbf{u} Y_k) = \nabla \cdot (D_k \nabla Y_k)
$$

When energy is enabled, the current passive sensible-enthalpy equation is:

$$
\frac{\partial h}{\partial t}
+
\nabla \cdot (\mathbf{u} h)
=
\frac{1}{\rho}\nabla \cdot (\lambda \nabla T)
+
\frac{q_{rad}}{\rho}
$$

where:

- `u` is velocity
- `p` is pressure
- `Y_k` are species mass fractions
- `rho` is constant density
- `nu` is kinematic viscosity (constant or Cantera-derived)
- `D_k` is species diffusivity (constant or Cantera-derived)
- `h` is transported sensible enthalpy
- `T = T(h,Y,p0)` is recovered temperature
- `lambda` is thermal conductivity
- `qrad` is a volumetric source term, currently zero unless a future coupling fills it
- `body_force` is used to drive periodic channel flow

---

## Mathematical Foundation: Projection Method

The solver uses a semi-implicit fractional-step method to decouple velocity and pressure:

1.  **Predictor Step**: Compute an intermediate velocity \(\mathbf{u}^*\) by advancing the momentum equation explicitly, excluding the new pressure gradient.
    $$\frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t} = -(\mathbf{u}^n \cdot \nabla) \mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n + \mathbf{f}$$
2.  **Poisson Equation**: Solve for the pressure potential \(\phi = p^{n+1} - p^n\) to enforce the divergence-free constraint.
    $$\nabla^2 \phi = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$$
3.  **Corrector Step**: Update the velocity and pressure fields.
    $$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla \phi$$
    $$p^{n+1} = p^n + \phi$$

---

## Current numerical method

The current flow method uses:

- cell-centered finite volume
- explicit momentum predictor
- AB2 explicit time integration after the first step
- forward Euler on the first step
- incremental pressure correction
- pressure Poisson solve using matrix-free CG/PCG
- conservative corrected face fluxes
- corrected cell-centered velocity update
- divergence diagnostic based on corrected face flux

Important recent numerical improvements:

- Conservative unnormalized pressure Poisson operator:
  - pressure matrix coefficient uses `area / distance`
  - pressure RHS uses `div * cell_volume`
  - this keeps the pressure matrix symmetric on nonuniform cuboid meshes

- Distance-weighted face interpolation:
  - replaces plain arithmetic averages for face values
  - reduces to the old `0.5 * (owner + neighbor)` form on uniform grids

- Periodic wrapped-distance handling:
  - periodic faces use paired-face geometry instead of raw long-distance cell-center separation
  - requires `periodic_face` and `periodic_neighbor` data from `periodic.dat`

- Pressure operator cache:
  - precomputes pressure neighbor IDs
  - precomputes pressure coefficients
  - precomputes pressure diagonal

- Persistent projection workspace:
  - avoids per-timestep allocation/deallocation of pressure and projection work arrays

Recommended flow settings for validation:

- use `convection_scheme = "central"` for laminar validation once stability is established
- use `convection_scheme = "upwind"` for first stability tests of new obstacle cases
- use periodic + body-force for channel flow
- use moving-wall + no-slip boundaries for lid-driven cavity
- avoid pressure inlet/outlet as a primary validation case until pressure boundary behavior is fully validated

---

## Current MPI architecture

Current MPI architecture:

- every rank reads and stores the full mesh
- full geometry is replicated on every rank
- flow computation is decomposed by contiguous owned cell ranges
- each rank computes only owned-cell contributions
- most global arrays are still replicated on every rank
- pressure matvec assembly uses owned-cell `MPI_Allgatherv` rather than full-array `MPI_Allreduce`
- rank 0 writes VTU/PVD/diagnostics output

Current decomposition:

```text
flow:
  decomposed by owned cell ranges

future radiation:
  decomposed by wavenumber / spectral interval

mesh:
  replicated on every rank
```

This architecture is intentional because the future radiation solver needs access to the full geometry on every rank.

Do not assume:

- distributed mesh
- graph partitioning
- PETSc distributed vectors
- parallel VTU/PVTU output
- fully local-only field storage

Current rule for matrix-free operators:

- geometry is globally available
- local computation fills owned entries only
- vectors used in neighbor access must be globally valid before the operator is applied
- pressure matvec output is assembled using owned-cell gather
- scalar dot products use global scalar reductions

---

## MPI and radiation design decision

The long-term radiation implementation will run every `N` flow steps and use the same MPI ranks as the flow solver.

Radiation will solve in the wavenumber/spectral domain:

- every rank has the full mesh
- every rank has the full scalar/species/radiation-relevant fields at radiation update time
- radiation work is decomposed by wavenumber
- each rank computes a partial `qrad(:)` contribution
- partial `qrad(:)` is reduced across ranks

Therefore, do not optimize the code by removing full mesh availability from ranks.

Preferred architecture:

```text
mesh:
  replicated globally

flow:
  owned-cell computation
  optimized replicated/global field communication for now
  possible future owned/ghost field layout

radiation:
  full geometry on each rank
  full radiation input fields on each rank
  wavenumber-domain decomposition
  MPI handled only in radiation driver
  radiation kernel is MPI-free
```

Future flow optimization should distinguish between:

```text
replicated mesh:
  keep

unnecessary full-vector communication every CG iteration:
  reduce or eliminate

full fields needed only every radiation step:
  gather/update only when needed
```

---

## Build commands

Debug build:

```bash
make clean
make BUILD=debug
```

Release build:

```bash
make clean
make BUILD=release
```

Recommended release flags:

```makefile
FFLAGS_RELEASE = -O3 -march=native -fno-omit-frame-pointer
```

Do not use `-ffast-math` during validation.

Aggressive performance flags may be tested only after validation:

```makefile
FFLAGS_RELEASE = -O3 -march=native -funroll-loops -fno-omit-frame-pointer -flto
```

Only keep aggressive flags if:

- diagnostics are unchanged
- velocity profiles are unchanged within tolerance
- serial and parallel trends still match
- runtime improves measurably

---

## Run commands

Run lid-driven cavity:

```bash
mpirun -np 1 ./lowmach_react_hex cases/lid_driven_cavity/case.nml
mpirun -np 8 ./lowmach_react_hex cases/lid_driven_cavity/case.nml
```

Run periodic channel:

```bash
mpirun -np 1 ./lowmach_react_hex cases/channel_flow/case.nml
mpirun -np 8 ./lowmach_react_hex cases/channel_flow/case.nml
```

If Makefile run targets support `NP`, prefer:

```bash
make cavity-release NP=8
make channel-release NP=8
```

Run square obstacle cases only after basic cavity/channel tests still pass.

---

## Required checks after solver changes

After modifying any of the following:

- flow solver
- pressure solver
- boundary conditions
- MPI routines
- mesh reader
- mesh generator
- output routines
- species routines
- energy/thermo routines
- radiation infrastructure
- input parsing
- mesh connectivity logic

Run these checks:

1. Build debug.
2. Build release.
3. Run cavity with 1 rank.
4. Run cavity with 8 ranks.
5. Run channel with 1 rank.
6. Run channel with 8 ranks.
7. If energy or thermo changed, run constant-property energy and Cantera thermo smoke tests.
8. Confirm:
   - no crash
   - VTU files are produced
   - PVD file points to the correct output files
   - `diagnostics.csv` is produced if enabled
   - `max_div` remains bounded
   - kinetic energy evolves smoothly
   - pressure iterations are not stuck at the maximum
   - serial and parallel trends match
   - output opens in ParaView

If the change affects pressure, also compare:

- `piter`
- pressure residual
- `max_div`
- `rms_div`
- `net_boundary_flux`
- kinetic energy

If the change affects mesh generation, also check:

- patch face counts
- periodic pair counts
- cell volume positivity
- face normal orientation
- `mesh%ncell_faces(c) == 6` for all current cuboid cells

---

## Current validation expectations

### Lid-driven cavity

For lid-driven cavity:

- Reynolds number is defined as `Re = U_lid * L / nu`.
- With `U_lid = 1`, `L = 1`, and `nu = 1e-2`, the case is `Re = 100`.
- `max_div` should remain bounded and small.
- Kinetic energy should evolve smoothly.
- For steady validation, use relative kinetic-energy change per output as a monitor.
- Compare centerline velocity profiles against benchmark data.

Steady-state indicator:

$$
\text{relative\_KE\_change} = \frac{|\text{KE}_{\text{new}} - \text{KE}_{\text{old}}|}{\max(|\text{KE}_{\text{new}}|, \text{tiny})}
$$

Typical target:

```text
relative_KE_change per output < 1e-5
```

Stricter validation target:

```text
relative_KE_change per output < 1e-6
```

The classic Ghia et al. benchmark is for a 2D cavity. For direct comparison, prefer a thin/symmetric/periodic spanwise direction rather than a fully 3D cavity with solid front/back walls.

### Periodic body-force channel

For periodic channel:

- x direction is periodic
- z direction is periodic
- y walls are no-slip
- flow is driven by `body_force_x`
- no pressure inlet/outlet is used

For a channel with walls at `y=0` and `y=H`, the steady analytic solution is:

$$
u(y) = \frac{\text{body\_force}_x}{2 \nu} y (H - y)
$$

Derived values:

$$
u_{\text{max}} = \frac{\text{body\_force}_x H^2}{8 \nu}
$$
$$
u_{\text{mean}} = \frac{\text{body\_force}_x H^2}{12 \nu}
$$

For channel Reynolds number, report explicitly:

```text
Re_H  = U_bulk * H / nu
Re_Dh = U_bulk * 2H / nu
```

For periodic body-force channel:

- `piter=0` can be acceptable when the predicted field is already divergence-free
- low body force should scale kinetic energy quadratically with velocity scale
- compare against the parabolic Poiseuille profile only after the flow reaches steady state

### Square/cube obstacle cases

For obstacle cases:

- first run low Reynolds number cases for stability
- use upwind initially
- switch to central only after the pressure solve is healthy
- monitor `piter`, `max_div`, `rms_div`, `net_boundary_flux`, and kinetic energy
- do not claim quantitative vortex-shedding accuracy without mesh/time refinement and benchmark comparison

Recommended square-obstacle validation ladder:

```text
Re = 20-40:
  steady symmetric wake

Re = 100:
  periodic shedding

Re = 250:
  stronger shedding after Re=100 is stable/validated
```

For square side length `D = 0.25` and `nu = 1e-2`:

```text
U(Re=40)  = 1.6
U(Re=100) = 4.0
U(Re=250) = 10.0
```

---

## Boundary-condition architecture

Current boundary-condition input supports:

- legacy `patch_type`
- separated `patch_velocity_type`
- separated `patch_pressure_type`

Do not remove legacy `patch_type`.

Example lid-driven cavity:

```fortran
&boundary_input
  n_patches = 6

  patch_name = "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"

  patch_type = "wall", "wall", "wall", "wall", "wall", "wall"

  patch_velocity_type = "no_slip", "no_slip", "no_slip", "moving_wall", "no_slip", "no_slip"

  patch_pressure_type = "zero_gradient", "zero_gradient", "zero_gradient", "zero_gradient", "zero_gradient", "zero_gradient"

  patch_u = 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
  patch_v = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  patch_w = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

  patch_p = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  patch_dpdn = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
/
```

Example periodic channel:

```fortran
&boundary_input
  n_patches = 6

  patch_name = "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"

  patch_type = "periodic", "periodic", "wall", "wall", "periodic", "periodic"

  patch_velocity_type = "periodic", "periodic", "no_slip", "no_slip", "periodic", "periodic"

  patch_pressure_type = "periodic", "periodic", "zero_gradient", "zero_gradient", "periodic", "periodic"

  patch_u = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  patch_v = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  patch_w = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

  patch_p = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  patch_dpdn = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
/
```

`dpdn` means:

```text
dpdn = dp/dn = grad(p) dot n
```

For zero-gradient pressure boundaries:

```text
patch_dpdn = 0.0
```

---

## Mesh architecture

Current solver mesh assumptions:

- hexahedral cells
- currently axis-aligned cuboid cells
- conformal mesh
- no hanging nodes
- no tetrahedra
- no prism/wedge cells
- no skewed/non-orthogonal hexahedra in the current FV operators

Native mesh files:

```text
points.dat
cells.dat
faces.dat
patches.dat
periodic.dat
```

Possible future explicit connectivity files:

```text
face_nodes.dat
cell_faces.dat
```

The current Fortran mesh reader reconstructs `mesh%ncell_faces` and `mesh%cell_faces` from `faces.dat`.

For current cuboid cells:

```text
mesh%ncell_faces(c) = 6
```

### Mesh generation strategy

Manual `.geo` files are acceptable for small tests but should not be the long-term workflow.

Recommended Level 1 mesh workflow:

```text
mesh_config.py
  -> native block mesh generator
  -> mesh_native/
  -> solver
```

or:

```text
mesh_config.json
  -> generated .geo
  -> gmsh
  -> convert_gmsh_hex.py
  -> mesh_native/
```

For the current solver, the direct native block generator is preferred for structured cuboid meshes.

The generator should:

- create axis-aligned cuboid cells
- assign owner/neighbor connectivity
- assign outward face normals
- assign patch IDs
- write periodic face pairing when requested
- verify face normal orientation
- verify every current cuboid cell has six faces
- print patch face counts

Keep Gmsh GUI as an inspection/debugging tool, not the only source of truth for production meshes.

---

## Current pressure-solver performance model

Current optimized pressure path:

- full mesh replicated
- global pressure vectors still available
- pressure coefficients cached
- pressure diagonal cached
- projection workspace cached
- pressure matvec fills owned entries
- owned entries are assembled with `MPI_Allgatherv`
- dot products use scalar `MPI_Allreduce`

This is more efficient than the earlier full-vector `MPI_Allreduce` after each pressure matvec.

Do not regress to full-array `MPI_Allreduce` for pressure matvec assembly unless explicitly debugging.

Near-term performance direction:

```text
keep full mesh replicated
keep radiation compatibility
reduce repeated allocation
cache mesh-dependent coefficients
avoid repeated MPI metadata setup
avoid full-vector reductions when owned-cell gather is sufficient
```

Long-term pressure options:

- stronger preconditioner
- geometric or algebraic multigrid
- PETSc/Hypre backend
- owned/ghost pressure vectors while keeping mesh replicated
- full distributed field architecture only after physics/validation stabilize

---

## Radiation architecture requirements

Radiation must be structured as:

```text
mod_radiation_api      # shared radiation data types
mod_radiation_driver   # MPI-aware wrapper
mod_radiation_kernel   # MPI-free implementation owned by radiation developer
```

The radiation developer should implement only:

```fortran
radiation_compute_qrad_kernel(mesh, inputs, spectrum, k_first, k_last, qrad)
```

The kernel must:

- not use `mpi_f08`
- not know rank count
- not perform MPI reductions
- receive full mesh
- receive full radiation-relevant scalar/species fields
- receive spectral data
- receive only assigned wavenumber range
- accumulate into `qrad(:)`

The driver handles:

- wavenumber partitioning
- MPI rank ownership of spectral intervals
- allocating partial `qrad`
- `MPI_Allreduce` or future `MPI_Reduce_scatterv`
- providing total `qrad` to flow/output

Radiation decomposition:

```text
mesh:
  replicated

fields needed by radiation:
  globally available at radiation update time

radiation work:
  decomposed by wavenumber

qrad:
  partial qrad per rank
  reduced to total qrad
```

Radiation call pattern:

```fortran
if (params%enable_radiation) then
   if (mod(step, params%radiation_interval) == 0) then
      call update_radiation_inputs(...)
      call radiation_update(...)
   end if
end if
```

Radiation should not affect the flow until the passive enthalpy path and manufactured `qrad` tests are validated. Until then, it may be computed and output for diagnostics/postprocessing.

---

## Cantera architecture requirements

Cantera should first be used as a property backend, not as a reaction source.

Do not call Cantera directly from:

- `mod_flow_projection`
- `mod_species_transport`
- tight flow/species loops
- radiation kernel

Preferred future modules:

- `mod_transport_properties`
- `mod_chemistry_cantera`
- `mod_thermo_state`

Suggested transport type:

```fortran
type transport_properties_t
   real(rk), allocatable :: rho(:)
   real(rk), allocatable :: mu(:)
   real(rk), allocatable :: nu(:)
   real(rk), allocatable :: lambda(:)
   real(rk), allocatable :: diffusivity(:,:) ! (nspecies,ncells)
end type transport_properties_t
```

Development path:

1. return constant transport properties (Done)
2. replace constants with Cantera-computed values (Done)
3. validate against known mixture values (Done)
4. Scale-on-Demand dynamic discovery (Done)
5. use transport properties in species and energy equations (Done)
6. use combined Cantera thermo sync for energy-step `T/cp/lambda/rho_thermo` updates (Done)

**Note**: Calling Cantera cell-by-cell at every step introduces significant overhead. The current energy thermo path uses a combined sync and conservative cache. Future work may add diagnostics for cache hit/miss behavior, tabulation, or stronger tolerance-based update strategies after validation.

---

## Stage 3A species requirements

Add passive transported species:

- CO2
- H2O
- CO

Initial model:

- no chemical reactions
- constant density
- constant viscosity
- constant species diffusivity
- species advected with corrected conservative `fields%face_flux`
- zero-gradient species at walls
- periodic species at periodic boundaries
- bounded mass fractions
- renormalize so `sum(Y_k)=1`
- output species to VTU:
  - `Y_CO2`
  - `Y_H2O`
  - `Y_CO`
  - `sum_Y`

Recommended species equation:

```text
dY_k/dt + div(u Y_k) = div(D_k grad(Y_k))
```

Use upwind advection for the first implementation to maintain boundedness. A limited second-order scheme can be added later.

Mixing is represented by the species diffusion term:

```text
div(D_k grad(Y_k))
```

For the first implementation:

- use constant species diffusivity
- later replace with Cantera-computed transport properties

Species validation tests:

1. Uniform species field should remain uniform.
2. Periodic scalar blob should advect without mass loss.
3. Species mass should be conserved for periodic/no-flux boundaries.
4. `sum(Y_k)` should remain near 1.
5. Serial and parallel species results should match.

---

## Stage 3E variable-density low-Mach formulation

Energy transport is now available as passive sensible enthalpy in the constant-density solver. The current projection constraint remains:

```text
div(u) = 0
```

Variable-density low-Mach reacting flow will generally require:

- density evolution
- equation of state
- energy/enthalpy equation
- thermodynamic pressure treatment
- modified divergence constraint
- coupling between heat release, species diffusion, radiation heat loss, and velocity divergence

This is a major architecture change and should come after:

- passive species works
- passive enthalpy works
- Cantera thermo sync is validated
- boundary enthalpy and species coupling tests pass
- chemistry source terms are understood
- radiation coupling through `qrad` is validated with manufactured sources

Do not patch variable-density behavior into the current incompressible projection without a written formulation.

---

## Output and postprocessing strategy

Short term:

- keep VTU/PVD for ParaView visualization
- keep diagnostics CSV
- add probe CSV files for quantitative analysis
- use Python scripts for validation profiles

Long term:

- VTU/PVD should remain the visualization/debugging format
- HDF5 should become the primary scientific data/restart/postprocessing format
- CSV should be used for reduced validation data

Recommended long-term output hierarchy:

```text
results.h5
  /mesh
    points
    cells
    cell_centers
    cell_volume
  /time
    step_000000
      time
      velocity
      pressure
      divergence
      Y_CO2
      Y_H2O
      Y_CO
      qrad
    step_000200
      ...
  /diagnostics
    step
    time
    kinetic_energy
    max_divergence
    rms_divergence
    pressure_iterations
```

Ideal workflow:

```text
solver writes:
  diagnostics.csv
  probes.csv
  results.h5
  results.xdmf
  optional flow_*.vtu for debugging/visualization

postprocessing reads:
  results.h5 for Python analysis
  results.xdmf or VTU/PVD for ParaView
```

Do not remove VTU output. It is still useful for debugging and ParaView.

---

## Coding style

- Use `implicit none`.
- Keep modules focused.
- Prefer clear names over clever names.
- Keep MPI-aware code isolated.
- Keep physics kernels MPI-free where practical.
- Prefer small subroutines with one responsibility.
- Avoid large rewrites unless the task explicitly asks for one.
- Use `fatal_error` for unrecoverable input/configuration errors.
- Keep rank-0-only file output unless parallel output is explicitly implemented.
- Keep public interfaces minimal.
- Avoid hidden changes to file formats.
- Update Makefile dependencies when adding modules.
- Preserve existing cases unless the task explicitly asks to change them.

---

## FORD documentation requirements

All new or significantly modified public modules, public types, and public procedures must have FORD-compatible documentation comments.

Use `!>` for documentation comments.

### Module documentation example

```fortran
!> MPI-aware radiation driver.
!!
!! This module owns all MPI logic for radiation. The radiation kernel must remain
!! MPI-free and should only compute local spectral contributions to qrad.
module mod_radiation_driver
```

### Public type documentation example

```fortran
!> Radiation input fields available on every rank.
!!
!! The arrays in this type are full-cell arrays, not owned-cell-only arrays,
!! because the radiation kernel needs full mesh visibility.
type, public :: radiation_inputs_t
   real(rk), allocatable :: temperature(:)  !< Cell temperature [K].
   real(rk), allocatable :: pressure(:)     !< Cell pressure [Pa].
end type radiation_inputs_t
```

### Public routine documentation example

```fortran
!> Compute total radiative source term qrad.
!!
!! The driver partitions the spectral grid by MPI rank, calls the MPI-free
!! radiation kernel for the local wavenumber range, and then reduces the partial
!! qrad contribution across ranks.
!!
!! @param mesh Full replicated mesh.
!! @param inputs Full replicated radiation input fields.
!! @param spectrum Spectral grid and quadrature data.
!! @param qrad Total cell-centered radiative source term.
subroutine radiation_update(mesh, inputs, spectrum, qrad)
```

### Documentation expectations

When adding a module, document:

- module purpose
- whether it is MPI-aware or MPI-free
- whether arrays are global, owned-only, or owned+ghost
- units of important physical quantities
- ownership of allocation/finalization
- validation assumptions

When adding a public routine, document:

- purpose
- input/output arguments
- MPI behavior
- assumptions
- whether it modifies global state
- whether it is safe for serial and parallel runs

When adding a derived type, document:

- what it stores
- array dimensions
- ownership semantics
- physical units where relevant

Codex should add comments that improve generated FORD documentation, not decorative comments that repeat obvious code.

### FORD / MathJax inline math delimiter rule

FORD/MathJax renders display equations correctly, but inline math in Markdown files and
Fortran documentation/source comments may be emitted literally when single-dollar
delimiters are used. Do not wrap inline math in single-dollar delimiters in documentation or source
comments. Use `\(...\)` for all inline math instead.

Examples of the required inline form:

```text
\(V_c\)
\(f \in \partial c\)
\(\rho_f\)
\(\psi_{\text{other}}\)
```

For Fortran source comments, write documentation like this:

```fortran
!> Cell control volume \(V_c\).
!! Face \(f \in \partial c\) uses density \(\rho_f\).
!! Other-side interpolation weight is \(\psi_{\text{other}}\).
```

Display equations may continue to use display-math delimiters where FORD/MathJax
renders them correctly. This rule is specifically about inline math. Before
committing documentation, grep Markdown and Fortran comments for single-dollar
inline math and replace those instances with `\(...\)`.

---

## Planned development sequence

Completed or active staged capabilities:

1. Passive species transport, no reactions.
2. Cantera transport-property interface.
3. Passive sensible-enthalpy energy transport.
4. Cantera thermo sync for `h <-> T`, `cp`, `lambda`, and diagnostic `rho_thermo`.
5. Profiling for projection, species, energy, Cantera sync, output, and MPI communication.

Near-term future work:

1. Manufactured `qrad` tests.
2. Radiation API/driver/kernel scaffold coupled through `qrad`.
3. Species-diffusion enthalpy correction.
4. Cantera chemistry source terms and heat release.
5. Variable-density low-Mach formulation.
6. Performance and scalable output.

