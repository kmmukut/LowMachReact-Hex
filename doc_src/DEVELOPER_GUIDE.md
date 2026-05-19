title: Developer Guide


# DEVELOPER_GUIDE.md — LowMachReact-Hex

This guide is for developers extending, validating, and maintaining
LowMachReact-Hex.  It describes the current codebase as it exists now, not only
the original constant-density baseline.

The solver should be understood as:

```text
stable baseline:
  transient constant-density incompressible / low-Mach finite-volume projection solver

active development layer:
  guarded non-reacting variable-density low-Mach solver using Cantera density,
  thermo, species, and transport coupling

future extensions:
  chemistry source terms, reaction heat release, radiation through qrad,
  restart/scientific output, and broader EOS/pressure validation
```

The solver is **not** a fully compressible acoustic Navier-Stokes code.  It is a
low-Mach/projection architecture where acoustic waves are filtered.

---

## 1. Project summary

LowMachReact-Hex is a Fortran 2008 MPI finite-volume solver for laminar flow on
hexahedral meshes.  The current code supports:

```text
- cell-centered finite-volume flow fields
- replicated global mesh on all MPI ranks
- owned-cell MPI decomposition
- fractional-step pressure projection
- passive species transport
- passive sensible-enthalpy transport
- Cantera-backed thermodynamics and transport
- active experimental non-reacting variable-density low-Mach mode
- VTU/PVTU/PVD visualization output
- CSV diagnostics and validation tooling
- profiling and Cantera cache statistics
```

Current development should preserve the stable constant-density path while
hardening the variable-density low-Mach path across more cases.

---

## 2. Current runtime modes

### 2.1 Constant-density baseline

Use this mode for hydrodynamic regression tests and first checks of new code.

```fortran
&fluid_input
  enable_variable_density = .false.
  density_eos = "constant"
  rho = ...
  nu = ...
/
```

Behavior:

```text
transport%rho = params%rho
projection target: div(u) = 0
mass flux: fields%mass_flux = params%rho * fields%face_flux
```

Cantera may still provide viscosity, diffusivity, cp, lambda, enthalpy, or
diagnostic thermo density, but flow density remains constant.

### 2.2 Cantera-assisted constant-density mode

This mode uses Cantera for properties while preserving the constant-density
projection.

Typical uses:

```text
- Cantera viscosity or species diffusivity
- Cantera h(T,Y,p0) and T(h,Y,p0)
- Cantera cp and lambda for energy transport
- diagnostic rho_thermo output
```

Important distinction:

```text
active flow density      = transport%rho
thermodynamic density    = energy%rho_thermo
```

When `enable_variable_density=.false.`, `rho_thermo` is diagnostic.

### 2.3 Experimental variable-density low-Mach mode

This mode is active but should remain under validation controls.

Expected configuration:

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

Behavior:

```text
transport%rho <- energy%rho_thermo
projection target: div(u) = S_projection
low-Mach source: S = (rho_old - rho)/(rho*dt) - (u.grad(rho))/rho
species update: conservative rho*Y branch
energy update: conservative rho*h branch
```

The mode is called experimental because validation is still being expanded
across boundary conditions, pressures, EOS choices, mesh/time resolution, and MPI
rank counts.

### 2.4 Not currently implemented

```text
- fully compressible acoustic Navier-Stokes
- shock-capturing compressible fluxes
- active chemical reaction source terms
- heat release from reactions
- coupled external radiation solver
- solver-side ideal-gas EOS path independent of Cantera
- full multicomponent Stefan-Maxwell FV diffusion operator
```

---

## 3. Repository structure

Typical developer-facing structure:

```text
src/
  Fortran solver modules and main program
  cantera_interface.cpp

tools/
  mesh conversion utilities
  diagnostics and validation scripts

doc_src/
  architecture, numerical method, input guide, validation, output, profiling docs

mechanisms/
  small Cantera thermo/transport demo mechanisms

cases/
  case setups, meshes, input files, and output directories
```

Generated output should normally stay under each case's configured output
directory:

```text
<output_dir>/VTK/
<output_dir>/diagnostics/
```

---

## 4. Major source modules

### 4.1 Main program

`src/main.f90` orchestrates:

```text
1. MPI startup
2. case namelist parsing
3. mesh import
4. transport initialization
5. flow/radiation MPI setup
6. boundary condition construction
7. field allocation and initialization
8. output preparation
9. transient timestep loop
10. finalization
```

The main loop is transient.  Steady solutions are obtained by time marching
until diagnostics stop changing.

### 4.2 Core infrastructure

| Module | Role |
|---|---|
| `mod_kinds` | precision, constants, fatal-error helpers |
| `mod_input` | namelist parsing, defaults, validation |
| `mod_mesh_types` | mesh/cell/face/patch data structures |
| `mod_mesh_io` | native mesh reader |
| `mod_fields` | flow field allocation and initialization |
| `mod_profiler` | hierarchical timing and profiling |

### 4.3 Flow and projection

| Module | Role |
|---|---|
| `mod_flow_projection` | momentum predictor, face fluxes, pressure solve, velocity/flux correction, CFL/diagnostics |
| `mod_bc` | field-specific boundary conditions |
| `mod_mpi_flow` | owned-cell decomposition and synchronization |

Projection targets:

```text
constant density:
  div(u) = 0

variable density:
  div(u) = S_projection
```

### 4.4 Species

`mod_species` owns passive species transport:

```text
- dynamic allocation for transported species
- advection by corrected face flux
- diffusion using constant or Cantera diffusivity
- correction velocity / sum(Y) control
- conservative rho*Y branch in variable-density mode
```

Reactions are not active in the current supported path.

### 4.5 Energy

`mod_energy` owns sensible enthalpy transport and thermodynamic dependent
fields:

```text
- h transport
- T recovery
- cp, lambda, rho_thermo storage
- qrad storage
- conduction
- species-enthalpy diffusion correction
- rho*h conservative branch in variable-density mode
- energy diagnostics and budget terms
```

The transported energy state is `h`, not `T`.

### 4.6 Transport and Cantera

| Component | Role |
|---|---|
| `mod_transport_properties` | constant/Cantera transport properties, active-density sync |
| `src/cantera_interface.cpp` | C ABI bridge to Cantera |

The Cantera bridge evaluates:

```text
mu, D_k
h_sens(T,Y,p0)
T(h,Y,p0)
cp
lambda
rho_thermo
species sensible enthalpies
```

`density_eos="cantera"` means active density comes from the selected Cantera
phase.  The actual EOS is defined inside the YAML phase, such as `ideal-gas`,
`Peng-Robinson`, or `Redlich-Kwong`.

### 4.7 Output and diagnostics

`mod_output` writes:

```text
- VTK visualization files
- flow diagnostics
- variable-density diagnostics
- compatibility diagnostics
- continuity residual diagnostics
- species/energy conservation diagnostics
- enthalpy budget diagnostics
- ParaView debug fields
```

CSV diagnostics are the source of truth for global validation.  VTK fields are
for spatial localization and qualitative debugging.

---

## 5. Mesh pipeline

The solver does not generate geometry internally.

```text
Gmsh .geo/.msh
  -> tools/mesh/convert_gmsh_hex.py
  -> points.dat / cells.dat / faces.dat / patches.dat / periodic.dat
  -> mod_mesh_io
  -> mesh_t
```

Current mesh assumptions:

```text
- hexahedral cells
- axis-aligned cuboidal volumes
- named physical boundary surfaces
- optional periodic face links through periodic.dat
```

Every MPI rank stores the full mesh.

---

## 6. MPI model

The current MPI model is intentionally simple:

```text
- full mesh replicated on every rank
- contiguous owned-cell ranges assigned by mod_mpi_flow
- owned cell updates only
- global synchronization through allreduce/allgather helpers
```

Benefits:

```text
- simple geometry access
- robust diagnostics
- easier radiation future path
- easier debugging and validation
```

Tradeoff:

```text
- replicated mesh limits very large mesh scalability
```

---

## 7. Boundary condition model

Boundary patches come from mesh patch names and are configured in `case.nml`.

Field-specific selectors:

```fortran
patch_type
patch_velocity_type
patch_pressure_type
patch_species_type
patch_temperature_type
```

Use field-specific selectors for new cases.  `patch_type` is a useful fallback
but should not be the only source of truth.

Accepted boundary concepts include aliases for:

```text
wall / no_slip / moving_wall
symmetry / symmetric / slip
periodic
dirichlet / fixed_value
neumann / zero_gradient
```

### 7.1 Current velocity-versus-mass-flux behavior

Current velocity boundaries prescribe velocity, which sets volumetric flux:

```text
prescribed u -> fields%face_flux
fields%mass_flux = rho_face * fields%face_flux
```

So in variable-density cases, equal velocities do not imply equal mass fluxes.

Future production-quality variable-density inlet handling should add:

```text
boundary density:
  rho_b = EOS(T_b,Y_b,p0)

mass-flow BC:
  target mdot -> u_n = mdot/(rho_b A)
```

### 7.2 Boundary thermodynamic state rule

For fixed-temperature and fixed-composition inlets:

```text
T_b = patch_T
Y_b = patch_Y
h_b = h(T_b,Y_b,p0)
rho_b = rho(T_b,Y_b,p0)
```

Do not use interior composition or a default bath composition when evaluating
fuel/oxidizer inlet thermodynamic states.

---

## 8. Cantera and EOS development rules

### 8.1 Mechanism and phase are case settings

Do not hard-code a mechanism or phase into solver logic or documentation.

```fortran
&fluid_input
  cantera_mech_file = "..."
  cantera_phase_name = "..."
/
```

Blank `cantera_phase_name` uses Cantera's default/first phase.  Nonblank selects
a named phase.

### 8.2 EOS is defined by the Cantera phase

The solver-level setting:

```fortran
density_eos = "cantera"
```

means density comes from Cantera.  The actual EOS is defined in the YAML phase:

```yaml
thermo: ideal-gas
thermo: Peng-Robinson
thermo: Redlich-Kwong
```

Transport model is also defined by the selected phase:

```yaml
transport: mixture-averaged
transport: multicomponent
transport: high-pressure
transport: high-pressure-Chung
```

### 8.3 Demo mechanisms are examples

The `mechanisms/` files are test/demo mechanisms for EOS and transport
selection.  They are not mandatory defaults.

Useful demo families:

```text
ideal_mixavg.yaml
ideal_multi.yaml
pr_high.yaml
pr_chung.yaml
rk_high.yaml
rk_chung.yaml
thermo_transport_demo.yaml
```

Use them to validate phase selection, property availability, and EOS/pressure
behavior before moving to detailed mechanisms.

### 8.4 Thermodynamic pressure is not projection pressure

Cantera uses:

```text
background_press = p0
```

The pressure projection solves a hydrodynamic pressure-like field.  Do not use
projection pressure boundary values as the Cantera thermodynamic pressure unless
the solver formulation is deliberately changed.

---

## 9. Energy and thermo rules

Core rule:

```text
h is transported
T is recovered from h,Y,p0
```

Sensible enthalpy convention:

```text
h_sens = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
```

Temperature recovery:

```text
h_abs_target = h_sens + h_abs(T_ref,Y,p0)
T = T(h_abs_target,Y,p0)
```

When species composition changes, preserve transported enthalpy:

```text
T_new = T(h_transported, Y_new, p0)
```

Do not rebuild enthalpy from old temperature:

```text
h_new = h(T_old, Y_new, p0)  ! wrong for the transported-state convention
```

The normal energy path should use the combined Cantera thermo sync when Cantera
thermo is enabled:

```text
(T, cp, lambda, rho_thermo) = sync(h,Y,p0)
```

---

## 10. Development workflow

### 10.1 Branching

Recommended branch policy:

```text
main:
  stable releases and production-ready code

dev:
  active integration branch

feature/<name>:
  isolated features or validation work

fix/<name>:
  targeted bug fixes
```

### 10.2 Patch discipline

For each numerical or diagnostic patch:

```text
1. State the purpose.
2. Identify whether it changes numerics, diagnostics, output, docs, or tooling.
3. Preserve constant-density regression behavior unless intentionally changed.
4. Add or update diagnostics when changing variable-density logic.
5. Update docs and case references.
6. Run the relevant validation subset.
```

### 10.3 Coding style

General expectations:

```text
- keep module responsibilities narrow
- avoid hiding numerical assumptions in utility routines
- prefer explicit mode checks over implicit side effects
- keep constant-density and variable-density branches readable
- preserve MPI ownership semantics
- use fatal_error or MPI-aware fatal paths for invalid solver states
- validate allocation sizes before use
```

### 10.4 When modifying Cantera logic

Before changing Cantera code:

```text
- identify whether the change affects constant-density or variable-density mode
- confirm T_ref convention
- confirm pressure is thermodynamic p0
- confirm composition source: cell Y, boundary Y, or bath/default
- update cache dependency rules if inputs change
- test mechanism/phase loading directly
```

---

## 11. Build and run workflow

Typical workflow:

```bash
# Load compiler/MPI/Cantera environment.

# Configure/build according to the project build system.
# Example only; use the repository's actual Makefile/CMake commands.
make clean
make

# Generate mesh if needed.
python tools/mesh/convert_gmsh_hex.py ...

# Run.
mpirun -np 4 ./lowmach_react_hex cases/<case>/case.nml
```

Output should appear under the case's configured `output_dir`.

Check:

```text
<output_dir>/VTK/
<output_dir>/diagnostics/
```

---

## 12. Output and postprocessing rules

### 12.1 Visualization

Use `.pvd` or `.pvtu` for full-domain ParaView visualization, especially for MPI
runs.  A single rank `.vtu` piece is not the full dataset.

Most solver arrays are cell-centered finite-volume fields.  ParaView filters
such as `Plot Over Line`, `Resample To Line`, or `Cell Data to Point Data` may
interpolate or smooth data.

Compare like with like:

```text
raw CellData extraction <-> raw FV cell-centered output
ParaView line sample    <-> equivalent sampled/interpolated script output
PVTU/PVD full dataset   <-> full-domain script extraction
```

### 12.2 Diagnostics

CSV diagnostics should be used for validation decisions:

```text
- projection residuals
- low-Mach source compatibility
- conservative continuity
- mass balance
- species integrals
- enthalpy/energy budget
- coupled transport audit
```

VTK/ParaView fields are for finding where a problem occurs.

---

## 13. Validation strategy

Validation should be a matrix, not one benchmark.

### 13.1 Always preserve baseline regressions

Before accepting changes, test at least one constant-density hydrodynamic case.
Variable-density changes should not break the stable baseline.

### 13.2 Flow validation

```text
- mesh import and volume/face consistency
- lid-driven cavity
- periodic body-force channel
- pressure outlet cases
- closed/wall-bounded pressure nullspace handling
- MPI rank-count consistency
```

### 13.3 Species validation

```text
- uniform species preservation
- scalar advection
- scalar diffusion
- boundedness and sum(Y)
- no-flux conservation
- fixed-composition inlet behavior
```

### 13.4 Energy validation

```text
- constant-property energy initialization
- Cantera T-h-T roundtrip
- conduction-only cases
- fixed-temperature boundary with boundary composition
- enthalpy advection
- species-enthalpy diffusion on/off
```

### 13.5 Variable-density validation

```text
- density sync from Cantera
- projection residual against S_projection
- conservative continuity residual
- mass-flux boundary diagnostics
- pressure/outlet compatibility
- variable_nu on/off
- dt refinement
- mesh refinement
- MPI rank-count changes
- EOS/pressure matrix
```

### 13.6 EOS and pressure validation

Include cases using:

```text
- ideal-gas Cantera phase
- Peng-Robinson Cantera phase
- Redlich-Kwong Cantera phase
- low background_press
- high background_press
- mixture-averaged transport
- high-pressure / Chung transport
```

### 13.7 Case-level comparisons

Counterflow is useful, but it is only one case family.  Other useful families:

```text
- closed-box density/temperature equilibration
- hot/cold channel
- rectangular mixing layer
- coflowing/opposed jets
- periodic scalar/thermal blobs
- manufactured qrad source tests
- high-pressure real-gas property boxes
```

---

## 14. Automated validation tools

Current validation tooling includes scripts for:

```text
- variable-density validation checks
- coupled transport conservation audit
- validation matrix runner
```

Use final-row CSV metrics for automated checks.  Primary metric categories:

```text
Projection:
  relative_divu_minus_S_projection_l2
  relative_divu_minus_S_projection_max

Conservative continuity:
  integral_drho_dt_plus_div_mass_flux_dV
  relative_conservative_residual_l2

Energy:
  relative_last_energy_update_balance_defect
  rel_output_recon_defect
  rel_operator_recon_defect
```

Warnings should be treated seriously even when a case does not hard-fail.

---

## 15. Profiling workflow

Enable profiling in the case file:

```fortran
&profiling_input
  enable_profiling = .true.
  nested_profiling = .true.
  write_cantera_cache_stats = .true.
/
```

Interpret nested timings carefully.  Inclusive rows are not additive when nested
timers are enabled.

Common profiling regions:

```text
Total_Simulation
Transport_Update
Projection_Step
Species_Transport
Energy_Transport
Diagnostics_Write_Flow
Diagnostics_Write_Energy
Output_Write_VTU
```

Cantera-heavy runs should also inspect cache statistics.

Do not optimize small timers before verifying correctness and conservation.

---

## 16. Common development mistakes

### 16.1 Treating variable-density mode as fully production-ready

The mode is active, but validation breadth is still growing.  Keep new changes
guarded and diagnostic-rich.

### 16.2 Confusing velocity BCs with mass-flow BCs

Velocity boundaries prescribe volume flux.  Mass flux is derived from density.
For controlled mass input, a mass-flow BC is a future separate feature.

### 16.3 Confusing thermodynamic pressure and projection pressure

`background_press` is for Cantera.  Projection pressure is for the pressure
solve.

### 16.4 Rebuilding enthalpy from old temperature after species transport

Preserve transported `h`; recover `T`.

### 16.5 Using the wrong composition at a boundary

Fixed-temperature species inlets need boundary `Y_b`, not interior `Y`.

### 16.6 Comparing ParaView interpolation to raw cell data

Use equivalent sampling/interpolation methods.

### 16.7 Reading one rank VTU from an MPI run

Use `.pvtu` or `.pvd` for full-domain analysis.

### 16.8 Assuming a Cantera phase's transport model equals FV operator physics

The FV species operator currently uses scalar diffusivities.  Full
multicomponent diffusion requires a different operator.

---

## 17. Near-term developer priorities

Recommended immediate work:

```text
1. Keep documentation synchronized with active variable-density behavior.
2. Add patch-wise mass-flow diagnostics for every physical boundary.
3. Add restart read/write support for long transient runs.
4. Centralize boundary thermodynamic state evaluation.
5. Expand the automated validation matrix across EOS, pressure, dt, mesh, and MPI.
6. Improve full-domain postprocessing scripts for PVTU/PVD output.
```

Short-to-medium term:

```text
1. Add density-aware boundary states.
2. Add optional mass-flow inlet BC.
3. Harden pressure/outlet compatibility in variable-density mode.
4. Add bounded higher-order scalar/energy convection option.
5. Add broader high-pressure EOS validation.
```

Longer term:

```text
1. Add non-reacting production validation coverage.
2. Add chemistry source terms and heat release.
3. Add manufactured and physical radiation coupling through qrad.
4. Consider HDF5/XDMF or similar for restart and scientific postprocessing.
5. Consider distributed mesh only if replicated mesh becomes limiting.
```

---

## 18. FORD source documentation

The project uses **FORD** for Fortran source documentation.  FORD should be
treated as part of the developer workflow, not as a separate afterthought.

### 18.1 What FORD documents

FORD is primarily used for the Fortran source tree:

```text
src/*.f90
```

It can extract:

```text
- module documentation
- public and private procedures
- derived types and fields
- call relationships where available
- comments written with FORD/Doxygen-style markup
```

The C++ Cantera bridge is documented separately in the markdown documentation
because FORD's main value in this project is Fortran module documentation.  The
bridge should still be linked from the generated documentation where practical.

### 18.2 How to write FORD-friendly comments

Use documentation comments on modules, derived types, public procedures, and
important internal helper routines.

Recommended style:

```fortran
!> One-line summary of the module or procedure.
!!
!! Longer description explaining the numerical or architectural role.
!! Include assumptions, mode-specific behavior, and conservation implications
!! when relevant.
module mod_example
```

For procedures:

```fortran
!> Advance one timestep of a transported scalar.
!!
!! @param mesh      Global mesh.
!! @param flow      Flow MPI ownership and communicator.
!! @param params    Runtime case parameters.
!! @param fields    Flow fields used for advection.
subroutine advance_example(mesh, flow, params, fields)
```

For derived types:

```fortran
!> Flow-field storage.
type :: flow_fields_t
   real(rk), allocatable :: u(:,:)      !< Cell-centered velocity [m/s].
   real(rk), allocatable :: p(:)        !< Projection pressure-like field.
end type flow_fields_t
```

### 18.3 Documentation expectations for new code

When adding or changing a module, update FORD-facing comments for:

```text
- module purpose
- public procedures
- derived types and allocatable fields
- mode-specific behavior: constant-density vs variable-density
- ownership assumptions: global mesh, owned cells, global arrays
- units and sign conventions
- diagnostic fields and CSV/VTK outputs
```

Do not document only the implementation mechanics.  Document the numerical
contract.  For example, a flux routine should state whether it returns
volumetric flux, mass flux, owner-outward flux, boundary-inward flux, or a
cell-centered divergence.

### 18.4 FORD and current-state documentation

FORD comments should describe the code as it currently behaves.  Avoid leaving
stale comments that say a feature is future-only after the code path becomes
active.

Particularly important areas to keep current:

```text
- variable-density active density sync
- low-Mach divergence source
- projection-source versus current-source time levels
- conservative rho*Y and rho*h branches
- Cantera density semantics
- boundary-condition assumptions
- output field meanings
```

The markdown documents should give the architectural explanation.  FORD comments
should make the source browsable and prevent developers from misusing routines.

### 18.5 Generating FORD documentation

The exact command depends on the repository's FORD project file.  Typical usage
is:

```bash
ford <ford-project-file>
```

or, if the repository uses a conventional project file name:

```bash
ford ford.md
```

Generated FORD output should not be committed unless the repository explicitly
tracks generated documentation.  Prefer committing the source comments and the
FORD configuration.


### 18.7 FORD math rules

Use math in FORD comments only when it improves the source documentation.  Keep
the math close to the code it documents and prefer short equations over long
derivations.

Recommended delimiters:

```text
Inline math:
  \( ... \)

Display math:
  \[
    ...
  \]
```

Examples:

```fortran
!> Computes the corrected low-Mach projection residual.
!!
!! The variable-density target is \( \nabla \cdot u = S \), so the
!! projection error is \( \nabla \cdot u - S \), not raw divergence.
```

```fortran
!> Low-Mach conservative density source.
!!
!! \[
!! S =
!! \frac{\rho^{old} - \rho}{\rho \Delta t}
!! -
!! \frac{u \cdot \nabla \rho}{\rho}
!! \]
```

Rules for maintainable FORD math:

```text
- Use \( ... \) for short inline symbols and equations.
- Use \[ ... \] for displayed equations.
- Keep equations ASCII in source comments where practical.
- Prefer \rho, \Delta t, \nabla, \cdot, \sum, and \partial over Unicode symbols.
- Do not use raw '<' or '>' in prose near formulas when Markdown may treat them
  as HTML; write "less than", "greater than", or use code formatting.
- Escape underscores in prose math names if they are not inside code spans or
  math delimiters.
- Put code identifiers in backticks, not math mode:
    `fields%mass_flux`, `transport%rho`, `energy%rho_thermo`
- Put physical variables in math mode:
    \( \rho \), \( u \), \( h \), \( Y_k \), \( p_0 \)
- Do not mix Fortran syntax and LaTeX syntax in the same expression.
- Avoid very long aligned derivations in source comments; place those in
  markdown method documentation instead.
```

Preferred notation for this project:

```text
Density:
  \( \rho \)

Velocity:
  \( u \)

Volumetric face flux:
  \( F_f = u_f \cdot n_f A_f \)

Mass face flux:
  \( \dot{m}_f = \rho_f F_f \)

Species:
  \( Y_k \)

Sensible enthalpy:
  \( h \)

Thermodynamic pressure:
  \( p_0 \)

Low-Mach source:
  \( S \)

Projection target:
  \( \nabla \cdot u = S \)

Conservative continuity:
  \( \partial_t \rho + \nabla \cdot (\rho u) = 0 \)
```

When documenting discrete finite-volume equations, clearly state the sign
convention.  For this code, face fluxes are generally owner-outward:

```fortran
!> Face mass flux is owner-outward.
!!
!! \[
!! \dot{m}_f = \rho_f (u_f \cdot n_f) A_f
!! \]
!!
!! A negative value on a physical boundary means inflow to the owner cell.
```

Do not rely on math notation alone for critical implementation details.  Always
pair equations with a short prose statement describing array names, units, and
time levels.


### 18.8 FORD review checklist

Before merging a change, check:

```text
- new public modules/procedures have FORD comments
- changed numerical behavior is reflected in comments
- units and sign conventions are documented
- variable-density behavior is not described as dormant if it is active
- Cantera bridge assumptions are linked or described in markdown docs
- generated FORD docs build without warnings that indicate broken references
```


---

## 19. Documentation maintenance rules

When changing code, update the relevant docs:

```text
architecture.md:
  high-level module/runtime architecture

current_solver_state.md:
  current physics and validation status

cantera_interface.md / cantera_cpp_source.md:
  Cantera bridge and interface behavior

input_configuration_guide.md:
  user-facing case.nml options and recommended configurations

numerical_method.md:
  governing equations and discretization

validation_metrics.md:
  accepted pass/fail metrics

paraview_output_fields.md:
  VTK/PVTU field meanings

output_layout.md:
  output directory contract

profiling.md:
  timers and performance interpretation
```

Avoid appending long historical patch logs into every document.  Keep current
docs concise, current-state oriented, and link to patch history only where
needed.

---

## 20. Status summary

```text
Constant-density flow solver:
  stable baseline

Cantera thermo/transport:
  active and modular

Passive species and sensible enthalpy:
  active

Variable-density non-reacting low-Mach mode:
  active, diagnostics-rich, still experimental pending broader validation

Reaction source terms:
  future work

Radiation coupling:
  future work through qrad

Restart:
  recommended near-term infrastructure addition

Fully compressible flow:
  outside current architecture
```