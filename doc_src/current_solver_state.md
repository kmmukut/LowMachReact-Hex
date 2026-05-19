
title: Current Solver State and Development Assessment


# Current Solver State and Development Assessment

This document describes the current solver state as a code-facing snapshot.  It
is intended to replace older notes that treated variable-density behavior as
future-only.  The current code has a stable constant-density path and
an active, guarded, non-reacting variable-density low-Mach path.  The latter is
still best treated as experimental until it is validated across a broader matrix
of cases, boundary conditions, EOS choices, pressures, mesh resolutions, and MPI
rank counts.

This document is deliberately not centered on one validation case.  Counterflow,
lid-driven cavity, periodic channel flow, scalar advection/diffusion, thermal
conduction, variable-density mixing, high-pressure EOS cases, and manufactured
source cases all serve different validation purposes.

---

## 1. One-line summary

The solver is a transient hexahedral finite-volume low-Mach/projection code with:

```text
- replicated global mesh on every MPI rank
- owned-cell flow decomposition
- field-specific boundary conditions
- constant-density incompressible baseline
- optional passive species transport
- optional passive sensible-enthalpy transport
- optional Cantera thermo/transport coupling
- guarded experimental non-reacting variable-density low-Mach mode
- diagnostics-heavy validation infrastructure
- scaffold for later radiation and chemistry coupling
```

It is **not** a fully compressible acoustic Navier-Stokes solver.  The current
direction is variable-density low-Mach reacting-flow capability, where acoustic
waves are filtered and the velocity field is constrained by a low-Mach
divergence condition.

---

## 2. Current runtime modes

### 2.1 Constant-density projection mode

This is the mature baseline.

Typical configuration:

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
face mass flux: fields%mass_flux = params%rho * fields%face_flux
Cantera rho_thermo, if computed, is diagnostic only
```

Cantera may still be used in this mode for viscosity, diffusivity, thermal
conductivity, heat capacity, sensible enthalpy, and temperature recovery, but
the active projection density remains the configured constant density.

### 2.2 Cantera-assisted constant-density mode

This mode keeps constant projection density while using Cantera for some or all
transport/thermo properties.

Useful combinations include:

```text
- constant rho and nu, Cantera species diffusivity
- constant rho, Cantera mu but fixed nu disabled or enabled by input
- constant rho, Cantera energy thermo sync T(h,Y,p0), cp, lambda, rho_thermo
```

The key distinction is:

```text
rho used by flow/projection = transport%rho
rho computed by Cantera      = energy%rho_thermo
```

When `enable_variable_density = .false.`, `energy%rho_thermo` should be treated
as an output/diagnostic thermodynamic density, not as the active continuity or
projection density.

### 2.3 Guarded experimental variable-density low-Mach mode

This mode is now active in the code and should no longer be documented as a
purely dormant scaffold.  It is guarded and non-reacting.

Expected configuration pattern:

```fortran
&fluid_input
  enable_variable_density = .true.
  density_eos = "cantera"
  enable_variable_nu = .true.   ! optional but common
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
active density source:
  transport%rho <- energy%rho_thermo

projection target:
  div(u) = S_projection

low-Mach source:
  S = (rho_old - rho)/(rho*dt) - (u.grad(rho))/rho

mass flux:
  fields%mass_flux = rho_f * fields%face_flux

species update:
  conservative rho*Y branch when variable-density mode is enabled

enthalpy update:
  conservative rho*h branch when variable-density mode is enabled
```

The mode should still be called **experimental** because its correctness depends
on continued validation of boundary conditions, pressure/outlet compatibility,
density/source time levels, scalar transport accuracy, EOS/pressure choices, and
MPI consistency.  However, it is no longer accurate to say that
`enable_variable_density=.true.` is only a parser scaffold.

---

## 3. Thermodynamics, EOS, and pressure model

### 3.1 Solver-level EOS selector

The solver-level selector is:

```fortran
density_eos = "constant" | "cantera" | "ideal_gas"
```

Current interpretation:

```text
constant:
  Use configured params%rho as active flow density.

cantera:
  Use the selected Cantera phase density as the active density when
  enable_variable_density = .true.; otherwise keep it diagnostic.

ideal_gas:
  Reserved for a future solver-side ideal-gas EOS path unless/until implemented.
```

### 3.2 Cantera phase selection

The solver should not hard-code a mechanism or phase.  A case chooses the
thermodynamic model through:

```fortran
&fluid_input
  cantera_mech_file = "path/to/mechanism.yaml"
  cantera_phase_name = "phase_name"   ! blank means Cantera default/first phase
/
```

The actual EOS and transport model are defined by the selected Cantera YAML
phase, for example:

```yaml
thermo: ideal-gas
transport: mixture-averaged
```

or:

```yaml
thermo: Peng-Robinson
transport: high-pressure
```

or:

```yaml
thermo: Redlich-Kwong
transport: high-pressure-Chung
```

For the current gas-flow path, treat gas-mixture phases as the validated target.
Cantera may support additional phase types, but surface, plasma, electrolyte,
condensed, and pure-fluid phases are not validated by this solver path unless a
specific case and validation plan are added.

### 3.3 Demo mechanisms are examples, not defaults

The `mechanisms/` directory contains curated demo mechanisms that can be selected
from a case file.  They are useful for controlled EOS/transport tests, but they
should not be described as required solver defaults.

Available one-phase demo files:

| File | Phase | Thermo/EOS model | Transport model |
|---|---:|---:|---:|
| `mechanisms/ideal_mixavg.yaml` | `gas` | `ideal-gas` | `mixture-averaged` |
| `mechanisms/ideal_multi.yaml` | `gas` | `ideal-gas` | `multicomponent` |
| `mechanisms/pr_high.yaml` | `gas` | `Peng-Robinson` | `high-pressure` |
| `mechanisms/pr_chung.yaml` | `gas` | `Peng-Robinson` | `high-pressure-Chung` |
| `mechanisms/rk_high.yaml` | `gas` | `Redlich-Kwong` | `high-pressure` |
| `mechanisms/rk_chung.yaml` | `gas` | `Redlich-Kwong` | `high-pressure-Chung` |

There is also a multi-phase demo file:

```text
mechanisms/thermo_transport_demo.yaml
```

with named phases such as `Ideal_MixAvg`, `PR_High`, and `RK_Chung`.  This file
is useful for checking the named-phase selector.

The YAML `state` block is only the phase's default load state.  During solver
execution, the code sets cell states from case-controlled temperature, pressure,
and mass fractions.

### 3.4 Thermodynamic pressure versus projection pressure

The solver uses two conceptually different pressures:

```text
background_press:
  thermodynamic pressure p0 passed to Cantera for h(T,Y,p0),
  T(h,Y,p0), rho(T,Y,p0), cp, lambda, and transport properties.

flow/projection pressure:
  hydrodynamic pressure-like field used by the projection method and by
  pressure boundary conditions.
```

Do not interpret fixed projection pressure at an outlet as direct replacement
for the thermodynamic pressure used by Cantera.  A high-pressure EOS case should
set the thermodynamic pressure through `background_press` and choose an
appropriate Cantera phase/mechanism.  Projection pressure boundary conditions
remain part of the low-Mach flow solve.

---

## 4. Main timestep ordering

The current high-level timestep logic is:

```text
0. Initialization
   - read case
   - read mesh
   - initialize transport
   - initialize MPI ownership
   - build boundary conditions
   - initialize flow/species/energy fields
   - if enabled, sync initial Cantera thermo density into active density

For each timestep:
1. Update CFL diagnostics / dynamic dt if enabled.
2. Refresh transport properties on transport_update_interval.
   - Cantera transport may use energy%T when energy is enabled.
   - In variable-density mode, transport refresh must not reset active rho
     back to params%rho.
3. Advance momentum and pressure projection.
   - constant density: div(u)=0
   - variable density: div(u)=S_projection
4. Advance species transport, if enabled.
5. Advance sensible enthalpy transport, if enabled.
   - preserve transported h
   - recover T(h,Y,p0)
   - refresh cp, lambda, rho_thermo
6. If variable-density mode is enabled:
   - sync transport%rho from energy%rho_thermo
   - update low-Mach source for the next projection
7. On output steps:
   - write flow diagnostics
   - write variable-density diagnostics if enabled
   - write energy/species/conservation diagnostics if enabled
   - write VTU/PVTU visualization files
```

This is a transient time-marching solver.  A steady solution is obtained only
when the time-marched fields stop changing according to chosen diagnostics.

---

## 5. Field and flux conventions

### 5.1 Volumetric and mass fluxes

The code distinguishes:

```text
fields%face_flux:
  corrected volumetric face flux, approximately u_f dot n_f A_f

fields%mass_flux:
  density-weighted face flux, approximately rho_f u_f dot n_f A_f
```

The boundary sign convention is outward-positive from the owner cell.  Inflow
through a physical boundary usually appears as negative outward flux.

### 5.2 Velocity boundary conditions prescribe volume flux

The current velocity boundary path prescribes boundary velocity, which sets a
volumetric face flux.  The solver then derives mass flux by multiplying by a
density:

```text
prescribed velocity -> face_flux
mass_flux = rho_face * face_flux
```

Therefore, in variable-density mode, equal inlet velocities generally do **not**
mean equal mass flux or equal momentum flux.  This is expected and applies to
all configurations, not only counterflow tests.

If a case requires a controlled mass inflow, a future density-aware mass-flow
boundary condition should set:

```text
target mdot -> face_flux = mdot / rho_boundary
```

rather than prescribing velocity directly.

### 5.3 Boundary density state

For physical boundary faces, current density weighting is primarily owner-cell
based.  For production-quality variable-density inlet/outlet treatment, boundary
density should eventually be computed from the boundary thermodynamic state:

```text
rho_b = EOS(T_b, Y_b, p0)
```

where `T_b` and `Y_b` come from the boundary condition when fixed-temperature or
fixed-composition inlets are used.

---

## 6. Boundary condition system

The solver supports field-specific boundary settings:

```fortran
patch_type
patch_velocity_type
patch_pressure_type
patch_species_type
patch_temperature_type
```

This split should be preserved.  It allows one patch to behave differently for
velocity, pressure, species, and temperature/enthalpy.

Common current interpretations:

```text
Velocity:
  no_slip / moving_wall / fixed_value / zero_gradient / symmetry / periodic

Pressure:
  zero_gradient / fixed_value / periodic / symmetry

Species:
  zero_gradient / fixed_value / periodic / symmetry

Temperature/enthalpy:
  fixed_value / zero_gradient / periodic / symmetry
```

The exact parser accepts aliases such as `dirichlet`, `fixed_value`,
`neumann`, `zero_gradient`, `wall`, `no_slip`, `symmetric`, and `periodic`.

Important boundary-state rule for energy:

```text
For fixed-temperature, fixed-composition inlets, boundary enthalpy must be
computed from h(T_b, Y_b, p0), not from the interior composition and not from a
default bath composition.
```

This rule is case-independent and should be part of the general validation plan.

---

## 7. Species transport state

The current species system supports passive, non-reacting mass-fraction
transport.  Reactions are intentionally disabled in the validated
variable-density path.

Current behavior:

```text
- species are advected by the corrected flow flux
- diffusion uses species diffusivities from constants or Cantera
- fields are bounded/renormalized so sum(Y_k) remains controlled
- variable-density mode uses a conservative rho*Y form
- correction-velocity style mass-conservation support is part of the design
```

Validation should include:

```text
- uniform species stays uniform
- scalar blob advection
- diffusion-only case
- periodic/no-flux species mass conservation
- fixed-composition inlet case
- serial vs MPI rank-count consistency
- species sum remains controlled
```

---

## 8. Energy and Cantera thermo state

The transported thermodynamic state is sensible enthalpy `h`, not temperature.

Core convention:

```text
h = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
T = T(h,Y,p0)
```

When species transport changes composition, the energy path should preserve the
transported enthalpy and recover the new temperature:

```text
T_new = T(h_transported, Y_new, p0)
```

It should not preserve old temperature by rebuilding:

```text
h_new = h(T_old, Y_new, p0)
```

Current energy capabilities:

```text
- sensible enthalpy transport
- Fourier conduction through grad(T)
- Cantera recovery of T from h,Y,p0
- Cantera cp, lambda, rho_thermo update
- optional species-enthalpy diffusion term
- qrad storage and diagnostics
- variable-density conservative rho*h update
```

Current missing physics:

```text
- chemical heat release
- reaction source terms
- external radiation model that fills qrad
```

---

## 9. Cantera transport limitations and interpretation

The Cantera bridge can load phases with different EOS and transport models.
However, the FV species operator currently consumes species diffusivity as
per-species scalar coefficients.  Even if a Cantera phase is configured with a
multicomponent transport model, the current solver does not yet solve a full
Stefan-Maxwell multicomponent diffusion system inside the FV operator.

Practical interpretation:

```text
mixture-averaged phase:
  natural match for the current scalar-D_k species diffusion operator.

multicomponent phase:
  useful for Cantera-supported properties and future development, but full
  multicomponent flux coupling requires additional bridge/operator work.

high-pressure or Chung transport:
  useful for real-gas/high-pressure testing if the selected Cantera phase can
  provide all properties the solver queries.
```

Before using a new mechanism or phase in production validation, probe it outside
the solver and verify at least:

```text
species names
thermo model
transport model
density
cp
thermal conductivity
viscosity
diffusivity availability
```

---

## 10. Variable-density diagnostics and validation metrics

Raw divergence is not the projection error in variable-density low-Mach mode.
The target is:

```text
div(u) = S_projection
```

Primary projection diagnostics should use:

```text
divu_minus_S_projection_max
divu_minus_S_projection_l2
relative_divu_minus_S_projection_max
relative_divu_minus_S_projection_l2
net_boundary_volume_flux_minus_integral_S_projection_dV
```

Current-source residuals are useful, but they measure source evolution after the
projection time level:

```text
divu_minus_S_current_*
S_current_minus_S_projection_*
net_boundary_volume_flux_minus_integral_S_current_dV
```

Primary conservative-continuity diagnostics should use:

```text
integral_drho_dt_plus_div_mass_flux_dV
conservative_residual_l2
relative_conservative_residual_l2
```

Primary energy/rho*h metrics should use the direct update closure and reconciled
budget metrics:

```text
relative_last_energy_update_balance_defect
rel_output_recon_defect
rel_operator_recon_defect
```

Older unreconciled or endpoint-only budget columns should be treated as
diagnostic context, not standalone failure metrics.

---

## 11. Output and postprocessing

The current output model is:

```text
<output_dir>/VTK/
  ParaView VTU/PVTU/PVD files

<output_dir>/diagnostics/
  CSV diagnostics and validation outputs
```

VTK data are mostly cell-centered finite-volume fields.  ParaView filters such
as `Plot Over Line`, `Resample To Line`, or `Cell Data to Point Data` may
interpolate or average fields.  Quantitative validation scripts should compare
like with like:

```text
raw FV CellData  <-> cell-center extraction
ParaView samples <-> equivalent interpolation in the script
PVTU/PVD dataset <-> full parallel dataset, not one rank piece
```

CSV diagnostics should be the source of truth for global conservation and
closure metrics.  VTK fields are best used for spatial localization and
qualitative debugging.

Long-term scientific output should move toward HDF5/XDMF or equivalent
restart/postprocessing data, while preserving VTU/PVTU for visualization.

---

## 12. General validation plan

Validation should be a matrix, not one case.

### 12.1 Baseline flow validation

```text
- mesh import and volume/face consistency
- serial vs MPI rank-count consistency
- lid-driven cavity
- periodic body-force channel
- pressure boundary condition tests
- periodic boundary tests
- moving-wall/no-slip/symmetry tests
```

### 12.2 Species validation

```text
- uniform species preservation
- periodic scalar advection
- diffusion-only species smoothing
- no-flux mass conservation
- fixed-composition inlet/outlet behavior
- boundedness and sum(Y_k)
```

### 12.3 Energy validation

```text
- Cantera T -> h -> T roundtrip
- fixed-temperature boundary with boundary composition
- conduction-only hot/cold wall case
- pure enthalpy advection
- species + enthalpy non-reacting coupling
- species-enthalpy diffusion on/off comparison
```

### 12.4 Variable-density low-Mach validation

```text
- density sync from Cantera thermo
- projection residual against S_projection
- conservative continuity residual
- mass-flux boundary balance
- variable_nu on/off
- ideal-gas vs real-gas Cantera phases
- low-pressure vs high-pressure background_press
- pressure boundary condition variants
- timestep refinement
- mesh refinement
- output-cadence variation
- MPI rank-count variation
```

### 12.5 Coupled transport validation

```text
- total mass trend
- transported species mass sum
- per-species boundary fluxes
- rho*h direct update closure
- reconciled output-state energy budget
- conduction boundary flux accounting
```

### 12.6 Case-level comparisons

Counterflow is one useful comparison because it exercises opposed inlets,
species gradients, enthalpy coupling, and density variation.  It should not be
the only validation target, and exact agreement with Cantera's 1D
`CounterflowDiffusionFlame` requires careful matching of model assumptions,
domain geometry, boundary conditions, mass fluxes, transport model, pressure,
and sampling procedure.

Other useful case families:

```text
- rectangular mixing layer
- opposed or coflowing jets
- closed-box variable-density diffusion
- hot/cold channel flow
- periodic scalar/thermal blobs
- manufactured qrad source cases
- high-pressure real-gas property boxes
```

---

## 13. Current code-level assessment

### 13.1 What is already strong

```text
- modular mesh/field/BC/projection/species/energy/transport separation
- replicated-mesh MPI model is simple and robust for current scale
- Cantera bridge isolates external thermodynamics from flow kernels
- named Cantera phase selector supports multiple EOS/transport choices
- variable-density diagnostics are unusually detailed
- energy budget diagnostics distinguish operator closure from bookkeeping
- output layout separates visualization from diagnostics
- validation tooling has started moving toward repeatable matrix checks
```

### 13.2 What needs hardening before calling variable-density production-ready

```text
- density-aware inlet/outlet boundary states
- optional mass-flow boundary condition or patch mass-flow controller
- clear compatibility rules for pressure outlets in variable-density mode
- broader EOS/pressure validation matrix
- grid/time refinement studies for scalar and energy transport
- MPI-rank-count regression across variable-density cases
- restart capability for long transient-to-steady runs
```

### 13.3 Numerical limitations to keep visible

```text
- first-order/upwind scalar transport can smear sharp species/temperature layers
- full multicomponent diffusion is not yet implemented in the FV operator
- boundary mass flux currently follows prescribed velocity and local density
- thermodynamic pressure is spatially uniform/background in the Cantera sync
- no reaction source terms or heat release are active
- no external radiation solver is coupled yet
- not a fully compressible acoustic solver
```

---

## 14. Recommended next development steps

### Near term

```text
1. Update stale documentation/comments to match the active variable-density path.
2. Add patch-wise mass-flow diagnostics for every physical boundary.
3. Add a small restart facility for long transient runs.
4. Add automated checks that compare PVTU/full-domain extraction, not single-rank VTU pieces.
5. Expand the validation matrix beyond the current baseline case.
```

### Short-to-medium term

```text
1. Add density-aware boundary-state evaluation:
   rho_b = EOS(T_b, Y_b, p0)

2. Add optional mass-flow inlet BC:
   target mdot -> u_n = mdot/(rho_b A)

3. Harden pressure/outlet compatibility for variable-density cases.

4. Improve scalar transport accuracy:
   bounded second-order / TVD / MUSCL option after first-order validation remains stable.

5. Extend Cantera transport support if true multicomponent FV diffusion is needed:
   pass and use coupled diffusion fluxes rather than only scalar D_k.
```

### Medium-to-long term

```text
1. Add non-reacting production validation across EOS, pressure, temperature,
   species contrast, mesh, dt, and MPI rank-count.

2. Add chemistry source terms and heat release only after non-reacting
   variable-density conservation is robust.

3. Couple radiation through qrad using manufactured qrad tests first.

4. Move restart/scientific postprocessing toward HDF5/XDMF while preserving
   VTU/PVTU visualization output.

5. Consider solver-side ideal-gas EOS only if it adds value beyond Cantera
   for performance, testing, or reduced-dependency runs.
```

---

## 15. Future goal assessment

Distance to major goals:

```text
Constant-density incompressible/low-Mach FV solver:
  implemented and should remain the stable baseline.

Non-reacting variable-density low-Mach solver:
  architecturally close and already active experimentally.
  Main remaining work is validation breadth, BC hardening, density-aware
  boundary states, restart, and scalar-transport accuracy.

Reacting variable-density low-Mach solver:
  requires chemistry source terms, heat release, reaction/transport
  time-integration choices, and stronger thermo/energy validation.

Radiating reacting low-Mach solver:
  requires external radiation physics coupled through qrad after manufactured
  qrad tests and energy-budget closure are stable.

Fully compressible Navier-Stokes solver:
  not the current architecture.  Would require a major formulation change,
  acoustic time integration or compressible fluxes, and a different pressure
  treatment.
```

The most realistic immediate target is:

```text
validated non-reacting variable-density low-Mach solver with Cantera EOS and
transport options, robust diagnostics, restart support, and density-aware
boundary conditions.
```

That target preserves the current architecture and gives a strong foundation for
later reactions and radiation.
