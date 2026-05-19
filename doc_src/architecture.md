title: Architecture


# LowMachReact-Hex Architecture

LowMachReact-Hex is a transient hexahedral finite-volume low-Mach/projection
solver intended as a clean foundation for non-reacting and eventually reacting
variable-density flow with Cantera thermodynamics/transport and future radiation
coupling.

The current architecture should be read as two supported layers:

```text
stable baseline:
  constant-density incompressible / low-Mach projection solver

active development layer:
  guarded non-reacting variable-density low-Mach solver using Cantera density,
  energy, species, and transport coupling
```

The code is not a fully compressible acoustic Navier-Stokes solver. Acoustic
compressibility is filtered by the projection/low-Mach formulation. The current
development target is a validated non-reacting variable-density low-Mach solver
that can later support reactions and radiation.

---

## 1. Design goals

The codebase is organized around these architectural goals:

```text
- keep geometry generation outside the solver core
- use simple, inspectable native mesh files
- keep full geometry available on every MPI rank
- separate flow, species, energy, transport, Cantera, output, and diagnostics
- support field-specific boundary conditions
- preserve a stable constant-density baseline
- add variable-density low-Mach capability incrementally and diagnostically
- keep thermodynamic EOS choice in the Cantera mechanism/phase when density_eos="cantera"
- provide validation-grade CSV diagnostics in addition to visualization output
- leave a clean qrad/radiation hook without entangling it with hydrodynamics
```

---

## 2. Mesh pipeline

The solver decouples mesh generation from the Fortran runtime.

```text
Gmsh geometry
  -> .msh file with named physical surfaces
  -> tools/mesh/convert_gmsh_hex.py
  -> native text mesh files
  -> mod_mesh_io
  -> mesh_t
```

### 2.1 Gmsh generation

Gmsh is used to generate structured or semi-structured hexahedral meshes with
named physical boundary surfaces. The solver relies on those physical names to
match case-file boundary condition entries.

### 2.2 Native conversion

`tools/mesh/convert_gmsh_hex.py` reads the `.msh` file, validates the supported
hexahedral topology, and writes the native mesh format:

```text
points.dat
cells.dat
faces.dat
patches.dat
periodic.dat
```

The current solver assumes axis-aligned cuboidal volume cells. The converter is
therefore part of the numerical contract, not just a convenience script.

### 2.3 Fortran import

`mod_mesh_io` reads the native files into `mesh_t`, whose major entities are:

```text
cells:
  volume, center, associated face ids

faces:
  owner, neighbor, periodic_neighbor, area, normal, center, patch id

patches:
  patch name, list of boundary face ids
```

The mesh is globally indexed and replicated on every MPI rank.

---

## 3. MPI ownership and decomposition

The solver uses a replicated-mesh, owned-cell MPI model.

### 3.1 Global mesh replication

Every MPI rank stores the full `mesh_t`. This is intentional:

```text
- it simplifies geometry access in matrix-free operators
- it simplifies diagnostics and output logic
- it supports future radiation work that may decompose by spectral/ray work
  rather than by spatial cells
```

This is not the most memory-scalable architecture for extremely large meshes,
but it is clear, robust, and appropriate for the current solver stage.

### 3.2 Flow-cell ownership

`mod_mpi_flow` assigns each rank a contiguous range of global cell ids:

```text
flow%first_cell : flow%last_cell
flow%owned(c)
```

Cell-centered updates are computed only for owned cells. Face-centered arrays
are stored globally, but face contributions are assembled consistently through
the owned-cell model and MPI synchronization routines.

### 3.3 Synchronization model

The flow communicator provides:

```text
global vector assembly:
  flow_allreduce_global_vector
  flow_allgather_owned_scalar

global reductions:
  MPI_Allreduce for dot products, norms, extrema, integrals, and diagnostics
```

Matrix-free pressure solves and global diagnostics use these synchronization
paths.

### 3.4 Separate radiation communicator

`mod_mpi_radiation` creates a separate communication layer for future radiation
tasks. The current hydrodynamic solver does not depend on a completed radiation
kernel, but the architecture keeps the door open for later work decomposition by
bands, rays, or external batches.

---

## 4. Runtime modes

### 4.1 Constant-density mode

Configuration pattern:

```fortran
&fluid_input
  enable_variable_density = .false.
  density_eos = "constant"
  rho = ...
/
```

Runtime behavior:

```text
active density:
  transport%rho = params%rho

projection target:
  div(u) = 0

mass flux:
  fields%mass_flux = params%rho * fields%face_flux
```

Cantera may still be used for viscosity, diffusivity, thermal conductivity,
heat capacity, enthalpy, and temperature recovery. In this mode, Cantera
thermodynamic density remains diagnostic and does not drive the pressure
projection.

### 4.2 Cantera-assisted constant-density mode

This mode keeps the constant-density projection while using Cantera for selected
thermo/transport properties.

Examples:

```text
- Cantera viscosity with or without variable_nu
- Cantera species diffusivities
- Cantera cp/lambda for energy transport
- Cantera h(T,Y,p0) and T(h,Y,p0)
- diagnostic rho_thermo output
```

The important separation is:

```text
flow density       = transport%rho
thermo density     = energy%rho_thermo
thermodynamic EOS  = selected Cantera phase
```

### 4.3 Experimental variable-density low-Mach mode

Configuration pattern:

```fortran
&fluid_input
  enable_variable_density = .true.
  density_eos = "cantera"
/

&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.
  thermo_update_interval = 1
/

&species_input
  enable_reactions = .false.
/
```

Runtime behavior:

```text
active density:
  transport%rho <- energy%rho_thermo

projection target:
  div(u) = S_projection

low-Mach source:
  S = (rho_old - rho)/(rho*dt) - (u.grad(rho))/rho

mass flux:
  fields%mass_flux = rho_f * fields%face_flux

species:
  conservative rho*Y branch

energy:
  conservative rho*h branch
```

This mode is active in the code but should still be called experimental until
the validation matrix covers more geometries, EOS choices, pressure levels,
boundary conditions, mesh/time refinements, and MPI rank counts.

### 4.4 Unsupported or not-yet-primary modes

The architecture currently does not provide:

```text
- fully compressible acoustic Navier-Stokes
- shock-capturing compressible fluxes
- active chemical source terms and heat release
- fully coupled radiation transport
- solver-side ideal-gas EOS path independent of Cantera
- production-grade full multicomponent Stefan-Maxwell diffusion inside the FV operator
```

---

## 5. Main program orchestration

`src/main.f90` owns the high-level simulation lifecycle:

```text
1. Start MPI.
2. Read case namelist.
3. Read native mesh.
4. Initialize transport fields.
5. Initialize flow and radiation MPI contexts.
6. Build boundary condition set.
7. Allocate and initialize flow fields.
8. Allocate and initialize species fields if enabled.
9. Allocate and initialize energy fields if enabled.
10. If variable-density mode is enabled, sync initial active density from thermo.
11. Prepare output directories and diagnostics.
12. Write initial output.
13. Enter transient timestep loop.
14. Finalize fields, mesh, communicators, and MPI.
```

The timestep loop is transient, not a direct steady-state nonlinear solve.
Steady solutions are obtained by running until field and diagnostic changes
become small.

---

## 6. Timestep data flow

The current timestep ordering is:

```text
A. CFL / dynamic timestep diagnostics
B. Cantera or constant transport-property update
C. Momentum prediction and pressure projection
D. Species transport
E. Sensible enthalpy transport
F. Cantera thermo sync
G. Variable-density active-density sync and low-Mach source update
H. Diagnostics and output on output steps
```

More explicitly:

```text
transport update:
  mu, nu, D_k, lambda as configured

projection:
  compute predicted velocity
  compute predicted face flux
  balance compatible outlet fluxes where appropriate
  solve pressure Poisson equation
  correct cell velocity and face flux
  compute mass flux

species:
  update Y or rho*Y depending on density mode

energy:
  update h or rho*h
  recover T from h,Y,p0
  refresh cp, lambda, rho_thermo

density:
  constant mode: preserve params%rho
  variable mode: transport%rho <- energy%rho_thermo

low-Mach source:
  assemble S for the next projection
```

The source used by a projection is saved separately as
`fields%projection_divergence_source`, because the current source may advance
after energy/thermo updates.

---

## 7. Major modules

### 7.1 Core infrastructure

| Module | Responsibility |
|---|---|
| `mod_kinds` | Precision, common constants, error handling helpers |
| `mod_input` | Namelist parsing, defaults, validation, runtime configuration |
| `mod_mesh_types` | Mesh, cell, face, and patch data structures |
| `mod_mesh_io` | Native mesh-file reader |
| `mod_fields` | Flow-field allocation and initialization |
| `mod_profiler` | Hierarchical MPI-aware timing diagnostics |

### 7.2 Flow and projection

| Module | Responsibility |
|---|---|
| `mod_flow_projection` | Momentum prediction, pressure equation, flux correction, CFL and flow diagnostics |
| `mod_bc` | Boundary condition mapping and field-specific boundary-state evaluation |
| `mod_mpi_flow` | Flow communicator, owned-cell decomposition, global assembly/reduction helpers |

The projection module supports the constant-density incompressible target and
the variable-density low-Mach target:

```text
constant density:
  div(u) = 0

variable density:
  div(u) = S_projection
```

The variable-density pressure operator uses density-dependent coefficients in
the active branch.

### 7.3 Species

| Module | Responsibility |
|---|---|
| `mod_species` | Passive species mass-fraction transport, diffusion, correction velocity, conservative variable-density branch |

Species are currently non-reacting. Reaction source terms are future work.

### 7.4 Energy

| Module | Responsibility |
|---|---|
| `mod_energy` | Sensible enthalpy transport, temperature recovery, cp/lambda/rho_thermo storage, qrad storage, energy diagnostics |

The transported energy state is sensible enthalpy:

```text
h = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
```

Cantera is used to recover:

```text
T = T(h,Y,p0)
cp = cp(T,Y,p0)
lambda = lambda(T,Y,p0)
rho_thermo = rho(T,Y,p0)
```

In variable-density mode, `rho_thermo` becomes the active flow density through
`mod_transport_properties`.

### 7.5 Transport and Cantera

| Module / file | Responsibility |
|---|---|
| `mod_transport_properties` | Constant or Cantera-backed transport properties; active-density sync |
| `src/cantera_interface.cpp` | C ABI bridge to Cantera C++ API |

The Cantera bridge handles:

```text
- mechanism and optional named-phase initialization
- species-name mapping
- mass-fraction vector construction
- viscosity and diffusivity evaluation
- sensible enthalpy evaluation
- temperature recovery from sensible enthalpy
- combined thermo sync for T, cp, lambda, rho_thermo
- cache statistics
```

`density_eos="cantera"` means that the active density is obtained from the
selected Cantera phase. The actual EOS is defined by that phase, not by the
string `"cantera"` itself.

### 7.6 Output and diagnostics

| Module | Responsibility |
|---|---|
| `mod_output` | VTU/PVTU/PVD visualization, CSV diagnostics, mesh summary, validation diagnostic files |

Output is separated as:

```text
<output_dir>/VTK/
  visualization files

<output_dir>/diagnostics/
  CSV diagnostics and validation files
```

Global conservation decisions should be based on diagnostics CSV files.
ParaView fields are primarily for spatial localization and qualitative analysis.

---

## 8. Boundary condition architecture

Boundary patches are defined in the mesh and configured in `case.nml`.

The system supports field-specific boundary selectors:

```fortran
patch_type
patch_velocity_type
patch_pressure_type
patch_species_type
patch_temperature_type
```

`patch_type` acts as a legacy or fallback selector. The field-specific entries
should be preferred for new cases.

### 8.1 Supported boundary concepts

The parser supports aliases corresponding to:

```text
wall / no_slip / moving_wall
symmetry / symmetric / slip
periodic
dirichlet / fixed_value
neumann / zero_gradient
```

These map into internal boundary condition ids used by the flow, pressure,
species, and temperature/enthalpy operators.

### 8.2 Field-specific interpretation

```text
velocity:
  fixed value, wall/no-slip, moving wall, symmetry/slip, zero-gradient,
  periodic

pressure:
  fixed value, zero-gradient, periodic, symmetry-compatible behavior

species:
  fixed composition, zero-gradient, periodic, symmetry-compatible behavior

temperature/enthalpy:
  fixed temperature, zero-gradient, periodic, symmetry-compatible behavior
```

### 8.3 Density-aware boundary-state gap

The current architecture primarily derives boundary mass flux from:

```text
prescribed boundary velocity -> volumetric face_flux
mass_flux = rho_face * face_flux
```

For variable-density production use, boundary density should eventually be
evaluated from the boundary thermodynamic state:

```text
rho_b = EOS(T_b, Y_b, p0)
```

This is especially important for fixed-temperature and fixed-composition inlets,
high-pressure EOS cases, and mass-flow-controlled boundaries.

### 8.4 Future mass-flow boundary condition

A mass-flow boundary condition would not replace velocity boundaries. It would
serve a different purpose:

```text
velocity BC:
  prescribe u, derive mdot = rho*u*A

mass-flow BC:
  prescribe mdot, derive u_n = mdot/(rho_b*A)
```

This is a natural future addition once boundary EOS state evaluation is
centralized.

---

## 9. EOS, pressure, and Cantera phase architecture

### 9.1 Solver-level density selector

The solver-level density selector is:

```text
density_eos = "constant" | "cantera" | "ideal_gas"
```

Current interpretation:

```text
constant:
  use params%rho

cantera:
  use the selected Cantera phase density when variable-density mode is enabled

ideal_gas:
  reserved for a future solver-side EOS implementation
```

### 9.2 Cantera mechanism and phase are case settings

The mechanism and phase are not architecture constants. They are case-level
configuration:

```fortran
&fluid_input
  cantera_mech_file = "..."
  cantera_phase_name = "..."
/
```

Blank phase name preserves Cantera's default/first phase behavior. A nonblank
phase name selects a specific phase from the YAML file.

The selected Cantera phase defines the actual thermodynamic and transport model:

```text
thermo: ideal-gas
thermo: Peng-Robinson
thermo: Redlich-Kwong

transport: mixture-averaged
transport: multicomponent
transport: high-pressure
transport: high-pressure-Chung
```

### 9.3 Thermodynamic pressure and projection pressure

The code distinguishes:

```text
background_press:
  thermodynamic pressure used for Cantera state calls

projection pressure:
  hydrodynamic pressure-like variable used by the pressure projection
```

High-pressure EOS tests should set `background_press` and choose an appropriate
Cantera phase. Projection pressure boundary conditions still control the
hydrodynamic pressure solve and should not be confused with the Cantera
thermodynamic pressure.

---

## 10. Transport architecture

### 10.1 Constant and Cantera-backed properties

Transport properties can come from configured constants or Cantera, depending
on input flags.

Tracked properties include:

```text
rho
rho_old
mu
nu
species diffusivity D_k
thermal conductivity lambda
cp
rho_thermo
```

In constant-density mode, active `rho` is preserved as the configured value. In
variable-density mode, active `rho` is synchronized from Cantera thermo after
energy sync.

### 10.2 Mixture-averaged versus multicomponent implications

The Cantera phase may expose mixture-averaged or multicomponent transport.
However, the current FV species operator consumes per-species diffusivity
coefficients. That naturally matches mixture-averaged transport.

A phase configured with multicomponent transport may be useful for future
development and some property calls, but a fully coupled multicomponent
Stefan-Maxwell FV diffusion operator is a separate architectural extension.

---

## 11. Energy and radiation architecture

### 11.1 Sensible enthalpy as transported state

The energy module transports sensible enthalpy `h`. Temperature is a dependent
thermodynamic variable recovered from:

```text
T = T(h,Y,p0)
```

This avoids rebuilding enthalpy from old temperature after species composition
changes.

### 11.2 Species-enthalpy diffusion

The code supports an optional species-enthalpy diffusion correction:

```text
-div(sum_k h_k J_k)
```

In variable-density mode, this path is density-weighted consistently with the
conservative `rho*h` update.

### 11.3 Radiation hook

`energy%qrad` is the current volumetric radiation/source hook. It is stored and
included in the energy equation/diagnostics, but no full external radiation
solver is currently coupled.

Future radiation coupling should:

```text
- fill qrad consistently at the correct time level
- validate manufactured qrad source cases first
- preserve energy-budget diagnostics
- use mod_mpi_radiation for non-spatial work decomposition when appropriate
```

---

## 12. Output and visualization architecture

### 12.1 VTU/PVTU/PVD

The solver writes XML VTK files for visualization:

```text
.vtu   per-rank piece files
.pvtu  parallel master file
.pvd   time-series collection
```

The PVTU schema must advertise every array that may appear in the per-rank VTU
pieces, including variable-density debug fields and energy reconciliation
fields.

### 12.2 Cell-centered finite-volume data

Most solver fields are finite-volume cell-centered fields. ParaView filters may
interpolate cell data to points or sample lines. Quantitative validation should
therefore compare equivalent representations:

```text
cell-centered extraction <-> raw CellData
line sampling            <-> same interpolation in analysis script
full PVTU/PVD dataset    <-> full-domain extraction, not one rank piece
```

### 12.3 CSV diagnostics

The diagnostics layer is a first-class architecture component. It provides:

```text
- flow diagnostics
- energy diagnostics
- variable-density projection diagnostics
- compatibility diagnostics
- continuity residual diagnostics
- species and enthalpy conservation diagnostics
- energy budget and reconciliation diagnostics
- coupled transport audit inputs
```

CSV diagnostics should remain the source of truth for conservation and closure
decisions.

---

## 13. Validation architecture

Validation is organized as a matrix, not a single benchmark.

### 13.1 Flow validation

```text
- mesh volume/face consistency
- pressure projection residuals
- closed-box no-flow preservation
- lid-driven cavity
- pressure outlet cases
- periodic channel/body-force flow
- MPI rank-count consistency
```

### 13.2 Species validation

```text
- uniform species preservation
- bounded scalar advection
- diffusion-only cases
- no-flux species mass conservation
- fixed-composition inlet behavior
- species sum control
```

### 13.3 Energy validation

```text
- Cantera T-h-T roundtrip
- conduction-only cases
- fixed-temperature boundary with boundary composition
- enthalpy advection
- species-enthalpy diffusion on/off
- qrad manufactured source cases
```

### 13.4 Variable-density validation

```text
- density sync from Cantera
- projection residual against S_projection
- conservative continuity residual
- mass-flux diagnostics
- pressure/outlet compatibility
- variable_nu on/off
- timestep refinement
- mesh refinement
- output-cadence variation
- MPI rank-count variation
```

### 13.5 EOS and pressure validation

```text
- ideal-gas Cantera phase
- Peng-Robinson Cantera phase
- Redlich-Kwong Cantera phase
- low-pressure and high-pressure background_press
- mixture-averaged and high-pressure transport phases
- mechanism/phase loading failures handled clearly
```

### 13.6 Case-level comparisons

Counterflow is one validation family, not the architecture center. Other useful
case families include:

```text
- closed-box density/temperature equilibration
- hot/cold channel flow
- periodic scalar and thermal blobs
- coflowing and opposed jets
- rectangular mixing layers
- pressure-driven variable-density channel cases
- manufactured-source qrad cases
- high-pressure real-gas property boxes
```

---

## 14. Known architectural limitations

The current architecture intentionally accepts some limitations:

```text
- globally replicated mesh limits very large mesh scalability
- axis-aligned cuboidal cell assumption limits geometry generality
- velocity BCs currently prescribe volume flux, not mass flux
- boundary density is not yet fully boundary-state/EOS based
- first-order/upwind scalar transport can smear sharp layers
- full multicomponent FV diffusion is not implemented
- reactions and heat release are not active
- radiation solver is not coupled
- restart capability is not yet a core runtime feature
- the solver is low-Mach/projection based, not fully compressible/acoustic
```

These limitations are acceptable for the current development target if they are
kept explicit in documentation and validation reports.

---

## 15. Future architecture direction

### 15.1 Near term

```text
- add patch-wise mass-flow diagnostics for all physical boundaries
- add simple restart read/write for long transient cases
- centralize boundary thermodynamic state evaluation
- expand automated validation matrices
```

### 15.2 Short-to-medium term

```text
- add density-aware inlet/outlet boundary states
- add optional mass-flow inlet boundary condition
- harden pressure/outlet compatibility for variable-density runs
- add bounded higher-order scalar/energy convection option
- add more EOS/pressure validation cases
- improve full-domain VTK/PVTU postprocessing tools
```

### 15.3 Medium-to-long term

```text
- add reaction source terms and heat release after non-reacting validation is robust
- add manufactured and physical radiation coupling through qrad
- consider HDF5/XDMF or similar for restart and scientific postprocessing
- consider distributed mesh storage only if mesh size requires it
- consider solver-side ideal-gas EOS only if it has a clear performance or
  dependency benefit over Cantera
```

---

## 16. Architectural status summary

```text
Constant-density low-Mach / incompressible solver:
  stable baseline

Cantera thermo/transport coupling:
  active and modular

Non-reacting variable-density low-Mach solver:
  active, diagnostics-rich, still experimental until broader validation hardens it

Reacting low-Mach solver:
  future extension after non-reacting variable-density conservation is robust

Radiating low-Mach solver:
  future extension through qrad and radiation communicator

Fully compressible solver:
  outside the current architecture
```
