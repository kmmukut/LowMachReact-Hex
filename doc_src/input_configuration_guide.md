title: Input Configuration Guide

# LowMachReact-Hex Input Configuration Guide

This guide explains how to configure LowMachReact-Hex simulations through the
Fortran namelist file `case.nml` and the native mesh files.

It is the long-form operating guide.  Use the compact `case_nml_reference.md`
when you only need a quick option table.

This version reflects the current implementation:

```text
stable baseline:
  transient constant-density incompressible / low-Mach projection solver

active guarded path:
  non-reacting variable-density low-Mach solver with Cantera thermo density,
  conservative species/enthalpy transport, and low-Mach source diagnostics

future extensions:
  reactions, heat release, full radiation coupling, restart/scientific output,
  and additional EOS/pressure validation
```

The word **guarded** means that the code path exists, but only for supported
input combinations.  It does not mean disabled.

---

## 1. Quick decision tree

### 1.1 Flow-only or constant-density validation

Use:

```fortran
enable_variable_density = .false.
density_eos = "constant"
```

Then:

```text
transport%rho = params%rho
projection target = div(u) = 0
rho_thermo, if present, is diagnostic
```

### 1.2 Cantera thermo/transport with constant projection density

Use this when you want Cantera properties but still want a constant-density flow
solve:

```fortran
enable_variable_density = .false.
density_eos = "constant"
enable_energy = .true.
enable_cantera_thermo = .true.
```

Then:

```text
flow/projection density = params%rho
thermodynamic density   = energy%rho_thermo, diagnostic only
T = T(h,Y,p0)
```

### 1.3 Non-reacting variable-density low-Mach mode

Use only the supported guarded combination:

```fortran
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
thermo_update_interval = 1
enable_reactions = .false.
```

Then:

```text
transport%rho <- energy%rho_thermo
projection target = div(u) = S_projection
species update    = conservative rho*Y branch
energy update     = conservative rho*h branch
```

`density_eos="cantera"` means "use density from the selected Cantera phase."
The actual EOS is the `thermo:` model inside the selected Cantera YAML phase.

### 1.4 Not yet supported as production modes

```text
enable_reactions = .true.              reaction source terms are future work
density_eos = "ideal_gas"              parsed/reserved, not the active solver-side EOS path
fully compressible acoustic flow        outside this solver architecture
full multicomponent FV diffusion         future operator extension
physical radiation solver               future; qrad hook exists
```

---

## 2. Recommended namelist block order

The solver reads named Fortran namelist blocks.  The recommended order is:

```fortran
&mesh_input
/

&time_input
/

&fluid_input
/

&solver_input
/

&boundary_input
/

&species_input
/

&energy_input
/

&output_input
/

&profiling_input
/
```

Namelist rules:

```text
- Strings must be quoted.
- Logical values are .true. and .false.
- Arrays may be comma-separated or explicitly indexed.
- Comments begin with !.
- Unknown names should be treated as input errors.
```

Example indexed boundary species input:

```fortran
patch_Y(1,1) = 0.20   ! species 1 on patch 1
patch_Y(2,1) = 0.05   ! species 2 on patch 1
patch_Y(3,1) = 0.75   ! species 3 on patch 1
```

---

## 3. Complete file skeleton

```fortran
&mesh_input
  mesh_dir = "cases/example/mesh_native"
/

&time_input
  nsteps = 1000
  dt = 1.0e-4
  output_interval = 100
  use_dynamic_dt = .false.
  max_cfl = 0.5
/

&fluid_input
  rho = 1.0
  nu = 1.0e-2

  enable_cantera = .false.
  enable_variable_density = .false.
  density_eos = "constant"
  enable_variable_nu = .false.

  cantera_mech_file = "gri30.yaml"
  cantera_phase_name = ""
  background_temp = 300.0
  background_press = 101325.0
  transport_update_interval = 10
/

&solver_input
  pressure_max_iter = 2000
  pressure_tol = 1.0e-8
  body_force_x = 0.0
  body_force_y = 0.0
  body_force_z = 0.0
  convection_scheme = "central"
/

&boundary_input
  n_patches = 0
/

&species_input
  enable_species = .false.
  enable_reactions = .false.
  enable_cantera = .false.

  nspecies = 0
  species_name = ""
  initial_Y = 1.0
  species_diffusivity = 0.0
/

&energy_input
  enable_energy = .false.
  enable_cantera_thermo = .false.
  enable_species_enthalpy_diffusion = .false.

  thermo_update_interval = 1
  thermo_default_species = "N2"

  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0

  energy_cp = 1005.0
  energy_lambda = 0.026
/

&output_input
  output_dir = "cases/example/output"
  write_vtu = .true.
  write_diagnostics = .true.
/

&profiling_input
  enable_profiling = .false.
  nested_profiling = .false.
  write_cantera_cache_stats = .false.
  variable_density_debug = .false.
/
```

`boundary_input` must be filled for real cases.  The empty block is only a
layout placeholder.

---

## 4. Native mesh input

`mesh_dir` points to native mesh files:

```text
points.dat
cells.dat
faces.dat
patches.dat
periodic.dat  only when periodic connectivity is used
```

### 4.1 `points.dat`

```text
npoints
id x y z
...
```

### 4.2 `cells.dat`

```text
ncells
id node1 node2 node3 node4 node5 node6 node7 node8 cx cy cz volume
...
```

### 4.3 `faces.dat`

```text
nfaces
id owner neighbor patch nx ny nz area cx cy cz
...
```

Meaning:

```text
owner     owner cell id
neighbor  adjacent cell id; boundary faces use a boundary/empty marker
patch     boundary patch id for physical boundary faces
nx ny nz  face normal components
area      face area
cx cy cz  face centroid
```

### 4.4 `patches.dat`

```text
npatches
id name nfaces
face_id_1 face_id_2 ...
...
```

Patch names must match `patch_name` in `case.nml`.

### 4.5 `periodic.dat`

Optional file for periodic connectivity:

```text
nlinks
face_id pair_face_id neighbor_cell_id
...
```

---

## 5. `&mesh_input`

| Option | Type | Required | Meaning |
|---|---:|---:|---|
| `mesh_dir` | string | yes | Directory containing native mesh files. |

Example:

```fortran
&mesh_input
  mesh_dir = "cases/channel_flow/mesh_native"
/
```

---

## 6. `&time_input`

| Option | Type | Required | Meaning |
|---|---:|---:|---|
| `nsteps` | integer | yes | Final timestep count. Must be non-negative. |
| `dt` | real | yes | Fixed timestep or initial timestep for dynamic CFL mode. Must be positive. |
| `output_interval` | integer | yes | Output/diagnostics cadence in timesteps. Must be positive. |
| `use_dynamic_dt` | logical | no | Enables CFL-based timestep adjustment. |
| `max_cfl` | real | if dynamic dt | Target CFL when dynamic timestep is enabled. |

Example:

```fortran
&time_input
  nsteps = 100000
  dt = 1.0e-4
  output_interval = 5000
  use_dynamic_dt = .false.
  max_cfl = 0.4
/
```

Interpretation:

```text
physical time advanced = nsteps * dt
```

when `use_dynamic_dt = .false.`.

---

## 7. `&fluid_input`

The `fluid_input` block controls active density, flow viscosity, Cantera phase
selection, thermodynamic pressure, and transport refresh.

```fortran
&fluid_input
  rho = 1.0
  nu = 1.0e-2

  enable_cantera = .false.
  enable_variable_density = .false.
  density_eos = "constant"
  enable_variable_nu = .false.

  cantera_mech_file = "gri30.yaml"
  cantera_phase_name = ""
  background_temp = 300.0
  background_press = 101325.0
  transport_update_interval = 10
/
```

### 7.1 `rho`

Constant/reference flow density `[kg/m^3]`.

Constant-density mode:

```text
transport%rho = params%rho
```

Variable-density mode:

```text
rho is an initialization/reference value
active density is later synced from energy%rho_thermo
```

For clean startup in variable-density cases, choose `rho` close to the expected
initial thermodynamic density.

### 7.2 `nu`

Configured kinematic viscosity `[m^2/s]`.

When variable viscosity is off, this is the flow viscosity used for momentum.

```fortran
enable_variable_nu = .false.
```

For fixed-Reynolds-number validation, keep `enable_variable_nu = .false.`.

### 7.3 `enable_cantera` in `fluid_input`

Controls the Cantera-backed fluid transport path, mainly viscosity.  It is
independent of `species_input enable_cantera`.

Use:

```fortran
&fluid_input
  enable_cantera = .true.
  enable_variable_nu = .true.
/
```

when you want Cantera viscosity to affect momentum.

Use:

```fortran
&fluid_input
  enable_cantera = .false.
  enable_variable_nu = .false.
/
```

when you want a fixed Reynolds number.

### 7.4 `enable_variable_density`

Activates the variable-density low-Mach path when used with a supported density
EOS and energy/thermo configuration.

Supported active configuration:

```fortran
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
thermo_update_interval = 1
enable_reactions = .false.
```

Do not enable variable density without Cantera thermo energy enabled.

### 7.5 `density_eos`

Solver-level active-density selector.

| Value | Current meaning |
|---|---|
| `"constant"` | Use configured `rho` as active flow density. |
| `"cantera"` | Use density from the selected Cantera phase when variable-density mode is enabled. |
| `"ideal_gas"` | Reserved for a future solver-side ideal-gas EOS path. |

Important:

```text
density_eos = "cantera"
```

does not mean "ideal gas".  It means the active density comes from Cantera.  The
actual EOS is defined by the selected YAML phase, such as `ideal-gas`,
`Peng-Robinson`, or `Redlich-Kwong`.

### 7.6 `enable_variable_nu`

Controls whether kinematic viscosity may vary with Cantera viscosity and active
density.

Typical choices:

```fortran
enable_variable_nu = .false.  ! fixed-Re validation
enable_variable_nu = .true.   ! variable-property flow validation
```

When `.true.`, fluid Cantera transport must be enabled.

### 7.7 `cantera_mech_file`

Mechanism path or Cantera-resolvable mechanism name.

Examples:

```fortran
cantera_mech_file = "gri30.yaml"
cantera_mech_file = "mechanisms/ideal_mixavg.yaml"
cantera_mech_file = "mechanisms/pr_high.yaml"
```

### 7.8 `cantera_phase_name`

Optional phase selector.

```fortran
cantera_phase_name = ""      ! Cantera default/first phase
cantera_phase_name = "gas"   ! named phase
```

Use a nonblank value when the mechanism contains multiple phases or when using
the demo one-phase files with phase name `gas`.

### 7.9 `background_temp`

Fallback transport temperature `[K]`.

```text
if energy is enabled:
  Cantera transport uses energy%T

if energy is disabled:
  Cantera transport uses background_temp
```

### 7.10 `background_press`

Thermodynamic pressure `[Pa]` passed to Cantera as `p0`.

Cantera uses this pressure for:

```text
h(T,Y,p0)
T(h,Y,p0)
cp(T,Y,p0)
lambda(T,Y,p0)
rho(T,Y,p0)
mu(T,Y,p0)
D_k(T,Y,p0)
```

This is not the projection pressure field.

### 7.11 `transport_update_interval`

Cadence for transport-property refresh.

Controls:

```text
fluid Cantera mu, if fluid Cantera path is enabled
species D_k, if species Cantera path is enabled
```

Does not control:

```text
Cantera h <-> T recovery
cp/lambda/rho_thermo energy sync
low-Mach source update
pressure projection
```

Recommended:

```text
1       validation/debugging
5-20    faster exploratory runs when property staleness is acceptable
```

---

## 8. `&solver_input`

```fortran
&solver_input
  pressure_max_iter = 2000
  pressure_tol = 1.0e-8
  body_force_x = 0.0
  body_force_y = 0.0
  body_force_z = 0.0
  convection_scheme = "central"
/
```

| Option | Type | Meaning |
|---|---:|---|
| `pressure_max_iter` | integer | Maximum pressure-solver iterations. |
| `pressure_tol` | real | Pressure-solver tolerance. |
| `body_force_x` | real | Body force in x direction. |
| `body_force_y` | real | Body force in y direction. |
| `body_force_z` | real | Body force in z direction. |
| `convection_scheme` | string | Momentum/scalar convection scheme selector. Common values: `"upwind"`, `"central"`. |

Recommendations:

```text
use "upwind" for first stability tests
use "central" for cleaner laminar validation once stable
```

For periodic channel flow, use body force rather than pressure inlet/outlet.

---

## 9. `&boundary_input`

Boundary input maps named mesh patches to field-specific boundary conditions.

```fortran
&boundary_input
  n_patches = 6
  patch_name = "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"

  patch_type = "dirichlet", "dirichlet", "neumann", "neumann", "symmetric", "symmetric"

  patch_velocity_type = "fixed_value", "fixed_value", "zero_gradient", "zero_gradient", "symmetric", "symmetric"
  patch_pressure_type = "zero_gradient", "zero_gradient", "fixed_value", "fixed_value", "symmetric", "symmetric"
  patch_temperature_type = "fixed_value", "fixed_value", "zero_gradient", "zero_gradient", "symmetric", "symmetric"
  patch_species_type = "fixed_value", "fixed_value", "zero_gradient", "zero_gradient", "symmetric", "symmetric"

  patch_u = 0.5, -0.5, 0.0, 0.0, 0.0, 0.0
  patch_v = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  patch_w = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

  patch_p = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  patch_dpdn = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

  patch_T = 300.0, 300.0, 300.0, 300.0, 300.0, 300.0

  patch_Y(1,1) = 0.20
  patch_Y(2,1) = 0.05
  patch_Y(3,1) = 0.75
/
```

### 9.1 Boundary options

| Option | Type | Meaning |
|---|---:|---|
| `n_patches` | integer | Number of configured patches. |
| `patch_name` | string array | Must match names in `patches.dat`. |
| `patch_type` | string array | Legacy/fallback patch type. Keep for compatibility. |
| `patch_velocity_type` | string array | Velocity boundary selector. |
| `patch_pressure_type` | string array | Pressure boundary selector. |
| `patch_temperature_type` | string array | Temperature/enthalpy boundary selector. |
| `patch_species_type` | string array | Species boundary selector. |
| `patch_u`, `patch_v`, `patch_w` | real arrays | Fixed or wall velocity components. |
| `patch_p` | real array | Fixed pressure value for pressure Dirichlet patches. |
| `patch_dpdn` | real array | Pressure normal gradient for Neumann patches. |
| `patch_T` | real array | Fixed temperature value for temperature Dirichlet patches. |
| `patch_Y(k,p)` | real array | Fixed species mass fraction for species `k` on patch `p`. |

### 9.2 Boundary type aliases

| Physical meaning | Accepted names |
|---|---|
| Wall | `"wall"`, `"no_slip"`, `"moving_wall"` |
| Symmetry/slip | `"symmetry"`, `"symmetric"`, `"slip"` |
| Periodic | `"periodic"` |
| Fixed value | `"dirichlet"`, `"fixed_value"` |
| Zero gradient | `"neumann"`, `"zero_gradient"` |

Use field-specific selectors for new cases.  Keep `patch_type` as a fallback.

### 9.3 Velocity boundaries

Current velocity boundaries prescribe velocity, therefore prescribe volumetric
face flux:

```text
given u_b -> face_flux = u_b dot n * area
mass_flux = rho_face * face_flux
```

In variable-density mode, equal inlet velocities do not imply equal mass fluxes.
A mass-flow boundary condition would be a separate future feature.

### 9.4 Pressure boundaries

Use `fixed_value` pressure only where needed to anchor or control the pressure
solve.  Use `zero_gradient` for walls/outlets when appropriate for the case.

For zero-gradient pressure:

```fortran
patch_pressure_type = "zero_gradient"
patch_dpdn = 0.0
```

### 9.5 Temperature boundaries

Setting `patch_T` alone is not enough.  For a fixed-temperature patch:

```fortran
patch_temperature_type = "fixed_value"
patch_T = ...
```

For a zero-gradient temperature patch, `patch_T` is not imposed.

### 9.6 Species boundaries

For fixed composition:

```fortran
patch_species_type = "fixed_value"
patch_Y(k,patch) = ...
```

Mass fractions on fixed-value species patches should sum to approximately 1.

### 9.7 Boundary thermodynamic state

For fixed-temperature and fixed-composition inlets:

```text
T_b = patch_T(p)
Y_b = patch_Y(:,p)
h_b = h(T_b,Y_b,p0)
rho_b = rho(T_b,Y_b,p0)
```

Dirichlet species boundaries use `patch_Y`.  Non-Dirichlet species boundaries
use the interior composition for the boundary state.

This matters for Cantera thermo, variable-density density, and future mass-flow
boundary support.

---

## 10. `&species_input`

```fortran
&species_input
  enable_species = .true.
  enable_reactions = .false.
  enable_cantera = .true.

  nspecies = 3
  species_name = "CH4", "O2", "N2"

  initial_Y = 0.10, 0.14, 0.76
  species_diffusivity = 2.0e-5, 2.0e-5, 2.0e-5
/
```

| Option | Type | Meaning |
|---|---:|---|
| `enable_species` | logical | Master switch for species transport. |
| `enable_reactions` | logical | Reaction source switch. Current supported validation path keeps this `.false.`. |
| `enable_cantera` | logical | Use Cantera species diffusivities `D_k`. |
| `nspecies` | integer | Number of transported species. |
| `species_name` | string array | Transported species names. Must match Cantera names when Cantera is used. |
| `initial_Y` | real array | Initial mass fractions. |
| `species_diffusivity` | real array | Constant fallback diffusivities `[m^2/s]`. |

### 10.1 `enable_reactions`

Current supported setting:

```fortran
enable_reactions = .false.
```

Do not use reactions until source terms and heat release are implemented and
validated.

### 10.2 `enable_cantera` in `species_input`

This is different from `fluid_input enable_cantera`.

```text
fluid_input enable_cantera:
  controls fluid transport, mainly viscosity

species_input enable_cantera:
  controls species diffusivity D_k
```

When `.false.`:

```text
D_k = species_diffusivity(k)
```

When `.true.`:

```text
D_k = Cantera diffusivity at T,Y,p0
```

Temperature source:

```text
energy enabled:  energy%T
energy disabled: background_temp
```

### 10.3 Species names

Species names must match the selected Cantera mechanism if any Cantera
thermo/transport path uses transported composition.

Example:

```fortran
species_name = "CH4", "O2", "N2"
```

---

## 11. `&energy_input`

```fortran
&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.
  enable_species_enthalpy_diffusion = .true.

  thermo_update_interval = 1
  thermo_default_species = "N2"

  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0

  energy_cp = 1005.0
  energy_lambda = 0.026
/
```

| Option | Type | Meaning |
|---|---:|---|
| `enable_energy` | logical | Master switch for sensible-enthalpy transport. |
| `enable_cantera_thermo` | logical | Use Cantera for `h <-> T`, `cp`, `lambda`, and `rho_thermo`. |
| `enable_species_enthalpy_diffusion` | logical | Include `-div(sum_k h_k J_k)` correction. Requires species + Cantera thermo for meaningful `h_k`. |
| `thermo_update_interval` | integer | Currently supported value is `1`. |
| `thermo_default_species` | string | Bath/default species when species transport is disabled. |
| `initial_T` | real | Initial temperature `[K]`. |
| `energy_reference_T` | real | Sensible enthalpy reference temperature `[K]`. |
| `energy_reference_h` | real | Constant-cp reference enthalpy `[J/kg]`. |
| `energy_cp` | real | Constant heat capacity fallback `[J/kg/K]`. |
| `energy_lambda` | real | Constant thermal conductivity fallback `[W/m/K]`. |

### 11.1 Energy state convention

The transported thermodynamic state is sensible enthalpy:

```text
h [J/kg]
```

Temperature is recovered:

```text
T = T(h,Y,p0)
```

When species changes composition, preserve `h` and recover the new temperature:

```text
T_new = T(h_transported,Y_new,p0)
```

Do not rebuild:

```text
h_new = h(T_old,Y_new,p0)
```

### 11.2 Constant-cp mode

When `enable_cantera_thermo = .false.`:

```text
h = energy_reference_h + energy_cp * (T - energy_reference_T)
T = energy_reference_T + (h - energy_reference_h) / energy_cp
lambda = energy_lambda
cp = energy_cp
```

### 11.3 Cantera thermo mode

When `enable_cantera_thermo = .true.`:

```text
h_sens = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
T = Cantera HP inverse using h_sens + h_abs(T_ref,Y,p0)
cp = cp(T,Y,p0)
lambda = thermal_conductivity(T,Y,p0)
rho_thermo = rho(T,Y,p0)
```

The preferred implementation path is the combined sync:

```text
(T, cp, lambda, rho_thermo) = sync(h,Y,p0)
```

The sync must preserve `h`.

### 11.4 `thermo_update_interval`

Current supported value:

```fortran
thermo_update_interval = 1
```

The thermodynamic state must remain logically synchronized every energy step.
Do not set this above 1 unless a stale-thermo strategy is implemented and
validated.

### 11.5 `thermo_default_species`

Used only when:

```text
enable_cantera_thermo = .true.
enable_species = .false.
```

Ignored when transported species are enabled.

Examples:

```fortran
thermo_default_species = "N2"
thermo_default_species = "CO2"
thermo_default_species = "H2O"
```

Use `"H2O"` with the letter `O`, not `"H20"` with zero.

### 11.6 Species-enthalpy diffusion

When enabled:

```text
-div(sum_k h_k J_k)
```

is included in the energy equation using species sensible enthalpies relative to
`energy_reference_T`.

Validate with both:

```fortran
enable_species_enthalpy_diffusion = .false.
enable_species_enthalpy_diffusion = .true.
```

---

## 12. `&output_input`

```fortran
&output_input
  output_dir = "cases/example/output"
  write_vtu = .true.
  write_diagnostics = .true.
/
```

| Option | Type | Meaning |
|---|---:|---|
| `output_dir` | string | Root output directory. |
| `write_vtu` | logical | Write VTU/PVTU/PVD visualization files. |
| `write_diagnostics` | logical | Write CSV diagnostics. |

Current output layout:

```text
<output_dir>/VTK/
  ParaView VTU/PVTU/PVD files

<output_dir>/diagnostics/
  CSV diagnostics and validation files
```

Before rerunning after changing output arrays, remove old output:

```bash
rm -rf cases/<case>/output
```

---

## 13. `&profiling_input`

```fortran
&profiling_input
  enable_profiling = .true.
  nested_profiling = .true.
  write_cantera_cache_stats = .true.
  variable_density_debug = .false.
/
```

| Option | Type | Meaning |
|---|---:|---|
| `enable_profiling` | logical | Enable MPI-aware profiler. |
| `nested_profiling` | logical | Show nested timer tree. Inclusive rows are not additive. |
| `write_cantera_cache_stats` | logical | Print final Cantera cache hit/miss statistics. |
| `variable_density_debug` | logical | Print verbose variable-density debug diagnostics. Use only for debugging. |

Common timer names:

```text
Transport_Update
Projection_Step
Species_Transport
Energy_Transport
Diagnostics_Write_Flow
Diagnostics_Write_Energy
Output_Write_VTU
```

Cantera-energy nested timers may include:

```text
Energy_Cantera_PreSync
Energy_PreFlux_Exchange
Energy_Flux_Update
Energy_Cantera_PostSync
Energy_Final_Exchange
```

---

## 14. Cantera mechanism and phase workflow

### 14.1 Demo mechanisms

Curated demo mechanisms may be used for EOS/transport tests:

| File | Phase | Thermo/EOS | Transport |
|---|---|---|---|
| `mechanisms/ideal_mixavg.yaml` | `gas` | `ideal-gas` | `mixture-averaged` |
| `mechanisms/ideal_multi.yaml` | `gas` | `ideal-gas` | `multicomponent` |
| `mechanisms/pr_high.yaml` | `gas` | `Peng-Robinson` | `high-pressure` |
| `mechanisms/pr_chung.yaml` | `gas` | `Peng-Robinson` | `high-pressure-Chung` |
| `mechanisms/rk_high.yaml` | `gas` | `Redlich-Kwong` | `high-pressure` |
| `mechanisms/rk_chung.yaml` | `gas` | `Redlich-Kwong` | `high-pressure-Chung` |

Example:

```fortran
&fluid_input
  density_eos = "cantera"
  cantera_mech_file = "mechanisms/pr_high.yaml"
  cantera_phase_name = "gas"
/
```

### 14.2 Multi-phase demo file

A multi-phase demo mechanism may contain phases such as:

```text
Ideal_MixAvg
Ideal_Multi
PR_High
PR_Chung
RK_High
RK_Chung
```

Select one with:

```fortran
cantera_phase_name = "PR_High"
```

Use exact case-sensitive phase names.

### 14.3 Probe a mechanism before using it

Run a quick Python probe before adding a mechanism to validation:

```bash
python - <<'PY'
import cantera as ct
gas = ct.Solution("mechanisms/pr_high.yaml", "gas")
gas.TPX = 300.0, ct.one_atm, {"CH4": 0.05, "O2": 0.21, "N2": 0.74}
print("thermo:", gas.thermo_model)
print("transport:", gas.transport_model)
print("rho:", gas.density)
print("cp:", gas.cp_mass)
print("lambda:", gas.thermal_conductivity)
print("mu:", gas.viscosity)
PY
```

### 14.4 Converting CHEMKIN files

Example:

```bash
ck2yaml --input chem.inp --thermo therm.dat --transport tran.dat --output mechanism.yaml
```

Then inspect the generated phase:

```bash
grep -n "^- name:" mechanism.yaml
grep -n "thermo:" mechanism.yaml | head
grep -n "transport:" mechanism.yaml | head
```

Use in `case.nml`:

```fortran
cantera_mech_file = "mechanism.yaml"
cantera_phase_name = "gas"
```

Detailed converted mechanisms are usually ideal-gas mechanisms unless manually
edited to include real-gas EOS data.  For `Peng-Robinson` or `Redlich-Kwong`,
every included species must have usable critical/EOS parameters.

---

## 15. Recommended configuration modes

### Mode 1: baseline hydrodynamics

```fortran
&fluid_input
  rho = 1.0
  nu = 1.0e-2
  enable_cantera = .false.
  enable_variable_density = .false.
  density_eos = "constant"
  enable_variable_nu = .false.
  background_temp = 300.0
  background_press = 101325.0
/

&species_input
  enable_species = .false.
  enable_reactions = .false.
  enable_cantera = .false.
  nspecies = 0
/

&energy_input
  enable_energy = .false.
  enable_cantera_thermo = .false.
/
```

Expected:

```text
constant density
div(u)=0
no species or energy fields required
```

### Mode 2: constant-property energy

```fortran
&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .false.
  enable_species_enthalpy_diffusion = .false.
  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/
```

Expected:

```text
h and T follow constant-cp relation
conduction uses grad(T)
rho remains constant
```

### Mode 3: Cantera thermo, no transported species

```fortran
&species_input
  enable_species = .false.
  nspecies = 0
/

&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.
  thermo_update_interval = 1
  thermo_default_species = "N2"
/
```

Expected:

```text
Y = bath/default species
T/cp/lambda/rho_thermo come from Cantera
rho_thermo is diagnostic unless variable-density mode is enabled
```

### Mode 4: Cantera thermo with transported species

```fortran
&species_input
  enable_species = .true.
  enable_reactions = .false.
  enable_cantera = .false.
  nspecies = 3
  species_name = "CH4", "O2", "N2"
  initial_Y = 0.10, 0.14, 0.76
/

&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.
  thermo_update_interval = 1
/
```

Expected:

```text
Cantera thermo uses transported Y
h is preserved when Y changes
T = T(h,Y,p0)
```

### Mode 5: Cantera species diffusivity with fixed-Re flow

```fortran
&fluid_input
  enable_cantera = .false.
  enable_variable_nu = .false.
  enable_variable_density = .false.
  density_eos = "constant"
/

&species_input
  enable_species = .true.
  enable_cantera = .true.
  enable_reactions = .false.
/
```

Expected:

```text
D_k = Cantera D_k(T,Y,p0)
nu = constant
rho = constant
Re = fixed
```

### Mode 6: guarded non-reacting variable-density low-Mach

```fortran
&fluid_input
  enable_variable_density = .true.
  density_eos = "cantera"
  enable_variable_nu = .true.
  enable_cantera = .true.
  cantera_mech_file = "mechanisms/ideal_mixavg.yaml"
  cantera_phase_name = "gas"
  background_press = 101325.0
/

&species_input
  enable_species = .true.
  enable_reactions = .false.
  enable_cantera = .true.
/

&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.
  thermo_update_interval = 1
/
```

Expected:

```text
transport%rho is synced from energy%rho_thermo
projection enforces div(u)=S_projection
rho*Y and rho*h conservative branches are active
reactions remain off
```

Use diagnostics to validate the run.  Do not judge variable-density projection
quality by raw `max_div`; use `divu_minus_S_projection_*`.

---

## 16. ParaView and output fields

### 16.1 Flow fields

Common fields:

```text
velocity
pressure
rho
nu
mass_flux_vector
mass_flux_divergence
```

`pressure` is projection pressure, not Cantera thermodynamic pressure.

### 16.2 Energy fields

Common fields:

```text
temperature
enthalpy
qrad
cp
thermal_conductivity
rho_thermo
thermal_diffusivity
thermo_pressure
```

Meaning:

| Field | Meaning |
|---|---|
| `rho` | Active flow density `transport%rho`. Constant in constant-density mode; Cantera-synced in variable-density mode. |
| `rho_thermo` | Thermodynamic density returned by Cantera. Diagnostic in constant-density mode; active density source in variable-density mode. |
| `thermo_pressure` | Uniform `background_press`. |
| `enthalpy` | Transported sensible enthalpy `h`. |
| `qrad` | Volumetric source. Positive adds energy to the gas. |

### 16.3 Species fields

Examples:

```text
Y_CH4
Y_O2
Y_N2
sum_Y
D_CH4
D_O2
D_N2
```

### 16.4 Variable-density debug fields

When variable-density mode is enabled, VTU/PVTU may include:

```text
lowmach_source_current
lowmach_source_projection
lowmach_source_difference
divu_recomputed
divu_minus_S_projection
divu_minus_S_current
rho_current
rho_projection
rho_current_minus_projection
mass_flux_divergence_recomputed
lowmach_source_history_estimate
lowmach_source_advective_density
u_dot_grad_rho
continuity_residual_estimate
rho_h_output_state
rho_h_operator_consistent
rho_h_density_reconciliation
relative_rho_h_density_reconciliation
```

Use these for spatial debugging.  Use CSV diagnostics for pass/fail decisions.

### 16.5 Cell data versus line sampling

Most fields are cell-centered finite-volume data.  ParaView filters such as
`Plot Over Line`, `Resample To Line`, and `Cell Data to Point Data` interpolate
or average fields.  Quantitative scripts should compare equivalent data:

```text
raw CellData  <-> raw cell-center extraction
line samples  <-> equivalent script interpolation
PVTU/PVD      <-> full-domain extraction, not one rank VTU piece
```

---

## 17. Diagnostics files

With diagnostics enabled, files are written under:

```text
<output_dir>/diagnostics/
```

Common files include:

```text
diagnostics.csv
energy_diagnostics.csv
variable_density_diagnostics.csv
variable_density_compatibility.csv
variable_density_transport_conservation.csv
variable_density_continuity_residual.csv
species_energy_conservation.csv
species_integrals.csv
enthalpy_energy_budget.csv
variable_density_worst_cell.csv
variable_density_worst_cell_faces.csv
variable_density_boundary_residual_summary.csv
variable_density_projection_audit_cells_rank<RANK>.csv
variable_density_projection_audit_faces_rank<RANK>.csv
```

Not every file appears in every mode.  Many are gated by:

```text
enable_variable_density = .true.
write_diagnostics = .true.
```

Primary variable-density projection metrics:

```text
divu_minus_S_projection_max
divu_minus_S_projection_l2
relative_divu_minus_S_projection_max
relative_divu_minus_S_projection_l2
```

Primary continuity metrics:

```text
integral_drho_dt_plus_div_mass_flux_dV
relative_conservative_residual_l2
```

Primary energy metrics:

```text
relative_last_energy_update_balance_defect
rel_output_recon_defect
rel_operator_recon_defect
```

---

## 18. Runtime reporting checklist

At startup, verify that the printed mode summary matches the intended run.

Check:

```text
Flow density mode
density_eos
Cantera mechanism file
Cantera phase name
background_press
flow viscosity mode
transport update interval
energy enabled/disabled
Cantera thermo enabled/disabled
species enabled/disabled
species Cantera enabled/disabled
variable-density mode enabled/disabled
variable-density debug enabled/disabled
```

For variable-density runs, also check that the runtime reporting makes clear:

```text
raw div(u) is not the projection error
the target is div(u)=S_projection
```

---

## 19. Validation ladder

### Step 1: flow-only regression

```fortran
enable_energy = .false.
enable_species = .false.
enable_variable_density = .false.
```

Expected:

```text
old hydrodynamic behavior unchanged
diagnostics.csv produced
VTK opens
```

### Step 2: species-only constant-density transport

```fortran
enable_species = .true.
enable_energy = .false.
enable_variable_density = .false.
```

Expected:

```text
sum_Y controlled
species remains bounded
species mass conserved in periodic/no-flux cases
```

### Step 3: constant-property energy

```fortran
enable_energy = .true.
enable_cantera_thermo = .false.
```

Expected:

```text
h/T constant-cp relation works
conduction-only cases behave smoothly
```

### Step 4: Cantera thermo without species

```fortran
enable_energy = .true.
enable_cantera_thermo = .true.
enable_species = .false.
```

Expected:

```text
thermo_default_species controls properties
T -> h -> T roundtrip is clean
```

### Step 5: Cantera thermo with species

```fortran
enable_species = .true.
enable_energy = .true.
enable_cantera_thermo = .true.
enable_reactions = .false.
```

Expected:

```text
Option A holds: preserve h, recover T from h,Y,p0
fixed-temperature species inlets use boundary Y
```

### Step 6: Cantera D_k with fixed-Re flow

```fortran
fluid_input enable_cantera = .false.
enable_variable_nu = .false.
species_input enable_cantera = .true.
```

Expected:

```text
D_k varies with T,Y,p0
flow Re remains fixed
```

### Step 7: variable-density non-reacting low-Mach

```fortran
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
enable_reactions = .false.
```

Expected:

```text
divu_minus_S_projection_* small
conservative continuity diagnostics acceptable
energy direct-update/reconciled metrics acceptable
```

### Step 8: manufactured qrad tests

Before a physical radiation model:

```text
uniform heating
uniform cooling
Gaussian source
localized source
```

Expected:

```text
qrad > 0 heats
qrad < 0 cools
integrated energy response matches integrated qrad
```

---

## 20. Common mistakes

### Mistake 1: Using stale variable-density documentation

The variable-density path is no longer just a future parser scaffold.  It is
active for the guarded non-reacting Cantera-density configuration.  Keep older
notes that say "`enable_variable_density=.true.` is unsupported" out of current
user-facing docs.

### Mistake 2: Enabling variable density without energy Cantera thermo

Wrong:

```fortran
enable_variable_density = .true.
enable_energy = .false.
```

Correct guarded pattern:

```fortran
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
thermo_update_interval = 1
enable_reactions = .false.
```

### Mistake 3: Confusing the two `enable_cantera` switches

```text
fluid_input enable_cantera:
  fluid viscosity path

species_input enable_cantera:
  species diffusivity path
```

### Mistake 4: Confusing projection pressure and thermodynamic pressure

```text
background_press = Cantera thermodynamic pressure p0
pressure field   = projection pressure-like variable
```

### Mistake 5: Setting `patch_T` without fixed temperature type

Wrong:

```fortran
patch_T = 400.0
```

Correct:

```fortran
patch_temperature_type = "fixed_value"
patch_T = 400.0
```

### Mistake 6: Species inlet mass fractions not summing to 1

Check:

```text
sum_k patch_Y(k,p) ~= 1
```

### Mistake 7: Using interior composition for inlet enthalpy

Correct for fixed-temperature fixed-composition inlets:

```text
h_b = h(T_b,Y_b,p0)
```

not:

```text
h_b = h(T_b,Y_interior,p0)
```

### Mistake 8: Recomputing enthalpy after species transport

Correct:

```text
preserve h
recover T = T(h,Y_new,p0)
```

Wrong:

```text
h = h(T_old,Y_new,p0)
```

### Mistake 9: Expecting equal velocity to mean equal mass flux

In variable-density mode:

```text
mdot = rho * u * A
```

Equal velocities with different densities produce different mass fluxes.

### Mistake 10: Judging variable-density projection by raw divergence

In variable-density low-Mach mode:

```text
div(u) = S
```

Use:

```text
divu_minus_S_projection_*
```

not raw `max_div`.

### Mistake 11: Reading one rank `.vtu` from an MPI run

For full-domain data, open:

```text
.pvtu
.pvd
```

not a single rank `.vtu` piece.

### Mistake 12: Reusing old output after schema changes

Remove output after adding/changing fields:

```bash
rm -rf cases/<case>/output
```

---

## 21. Minimal examples

### 21.1 Minimal hydrodynamic case

```fortran
&mesh_input
  mesh_dir = "cases/channel_flow/mesh_native"
/

&time_input
  nsteps = 1000
  dt = 1.0e-4
  output_interval = 100
  use_dynamic_dt = .false.
  max_cfl = 0.5
/

&fluid_input
  rho = 1.0
  nu = 1.0e-2
  enable_cantera = .false.
  enable_variable_density = .false.
  density_eos = "constant"
  enable_variable_nu = .false.
  cantera_mech_file = "gri30.yaml"
  cantera_phase_name = ""
  background_temp = 300.0
  background_press = 101325.0
  transport_update_interval = 10
/

&solver_input
  pressure_max_iter = 2000
  pressure_tol = 1.0e-8
  body_force_x = 0.0
  body_force_y = 0.0
  body_force_z = 0.0
  convection_scheme = "central"
/

&boundary_input
  n_patches = 0
/

&species_input
  enable_species = .false.
  enable_reactions = .false.
  enable_cantera = .false.
  nspecies = 0
/

&energy_input
  enable_energy = .false.
  enable_cantera_thermo = .false.
  enable_species_enthalpy_diffusion = .false.
  thermo_update_interval = 1
  thermo_default_species = "N2"
  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/

&output_input
  output_dir = "cases/channel_flow/output"
  write_vtu = .true.
  write_diagnostics = .true.
/

&profiling_input
  enable_profiling = .false.
  nested_profiling = .false.
  write_cantera_cache_stats = .false.
  variable_density_debug = .false.
/
```

### 21.2 Guarded variable-density non-reacting case pattern

```fortran
&fluid_input
  rho = 1.0
  nu = 1.5e-5
  enable_cantera = .true.
  enable_variable_density = .true.
  density_eos = "cantera"
  enable_variable_nu = .true.
  cantera_mech_file = "mechanisms/ideal_mixavg.yaml"
  cantera_phase_name = "gas"
  background_temp = 300.0
  background_press = 101325.0
  transport_update_interval = 1
/

&species_input
  enable_species = .true.
  enable_reactions = .false.
  enable_cantera = .true.
  nspecies = 3
  species_name = "CH4", "O2", "N2"
  initial_Y = 0.10, 0.14, 0.76
  species_diffusivity = 2.0e-5, 2.0e-5, 2.0e-5
/

&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.
  enable_species_enthalpy_diffusion = .true.
  thermo_update_interval = 1
  thermo_default_species = "N2"
  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/
```

---

## 22. Documentation maintenance note

Do not append long historical patch logs to this guide.  Keep this document as a
current operating guide.

Patch history belongs in a changelog or development notes.  User-facing input
documentation should answer:

```text
What parameter exists?
What values are valid?
What does it do now?
What combinations are supported?
What diagnostics confirm it worked?
```