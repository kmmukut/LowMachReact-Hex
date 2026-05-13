title: Input Configuration Guide

# LowMachReact-Hex Input Configuration Guide

This guide explains how to configure LowMachReact-Hex simulations using the `case.nml` Fortran namelist file and native mesh input files.

It updates the earlier draft for the current energy/species/Cantera stage. The companion file `case_nml_reference.md` is the concise option reference. This file is the longer operating guide with recommended modes, examples, validation steps, and common mistakes.

## 1. Current development assumptions

The current solver should be understood as a constant-density low-Mach/incompressible projection solver with optional species and enthalpy transport.

Important current rules:

1. The pressure projection uses constant flow density `rho`.
2. Cantera thermodynamic density `rho_thermo` is diagnostic/future-use only.
3. `enable_variable_density = .true.` is not supported yet.
4. Energy transport uses sensible enthalpy `h` as the transported state.
5. Temperature is recovered from `h`, composition `Y`, and thermodynamic pressure `p0`.
6. `p0` is currently the uniform operating pressure `background_press`.
7. Cantera thermodynamics use sensible enthalpy, not raw absolute enthalpy.
8. Option A is the accepted species-energy coupling convention: after species changes composition, preserve transported enthalpy and recover the new temperature.
9. Fixed-temperature species inlets should use boundary composition when evaluating boundary enthalpy.
10. Radiation physics is not coupled yet, but `qrad` exists as the future volumetric source interface.

## 2. How the solver reads `case.nml`

The solver reads named Fortran namelist blocks from the case file. Recommended block order is:

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

The order is chosen to match the logical setup sequence: mesh, time, fluid properties, solver controls, boundary conditions, species, energy, output, and profiling.

General Fortran namelist rules:

- Strings should be quoted.
- Logical values are `.true.` and `.false.`.
- Arrays can be specified as comma-separated lists or with explicit indices.
- Comments begin with `!`.
- Unknown names in a namelist block should be treated as input errors.

Example indexed species boundary input:

```fortran
patch_Y(1,1) = 0.0   ! O2 on patch 1
patch_Y(2,1) = 0.0   ! N2 on patch 1
patch_Y(3,1) = 1.0   ! CO2 on patch 1
```

## 3. Native mesh input files

`mesh_dir` points to a directory containing the solver native mesh files.

### `points.dat`

```text
npoints
id x y z
...
```

### `cells.dat`

```text
ncells
id node1 node2 node3 node4 node5 node6 node7 node8 cx cy cz volume
...
```

### `faces.dat`

```text
nfaces
id owner neighbor patch nx ny nz area cx cy cz
...
```

Meaning:

- `owner`: owner cell id.
- `neighbor`: neighboring cell id. Boundary faces usually use a boundary/empty neighbor indicator.
- `patch`: boundary patch id for boundary faces.
- `nx ny nz`: face normal components.
- `area`: face area.
- `cx cy cz`: face centroid.

### `patches.dat`

```text
npatches
id name nfaces
face_id_1 face_id_2 ...
...
```

Patch names in this file must match `patch_name` in `case.nml`.

### `periodic.dat`

Optional file used only for periodic connectivity.

```text
nlinks
face_id pair_face_id neighbor_cell_id
...
```

## 4. Fluid configuration

The `&fluid_input` block controls the flow density, flow viscosity, thermodynamic operating pressure, and Cantera transport refresh interval.

```fortran
&fluid_input
  rho = 1.0
  nu = 1.0e-2

  enable_cantera = .false.
  enable_variable_density = .false.
  enable_variable_nu = .false.

  cantera_mech_file = "gri30.yaml"
  background_temp = 300.0
  background_press = 101325.0

  transport_update_interval = 10
/
```

### `rho`: constant flow density

`rho` is the density used by the flow/projection solver and by the constant-density energy equation.

Current behavior:

```text
flow density = params%rho
```

This is not replaced by Cantera density. The Cantera density field is diagnostic only.

### `nu`: constant kinematic viscosity

`nu` is the constant kinematic viscosity used when variable Cantera viscosity is not active.

For fixed-Reynolds-number validation, keep:

```fortran
enable_variable_nu = .false.
```

Then the Reynolds number is controlled by:

```text
Re = U_ref L_ref / nu
```

### `enable_cantera` in `fluid_input`

This switch controls the Cantera fluid transport path, mainly viscosity. It is independent of `species_input enable_cantera`.

Recommended current validation setting:

```fortran
enable_cantera = .false.
enable_variable_nu = .false.
```

Turn it on only for a dedicated variable-viscosity validation.

### `enable_variable_density`

This is a future variable-density low-Mach flag. Current supported value:

```fortran
enable_variable_density = .false.
```

Do not enable it yet. Cantera `rho_thermo` is available for diagnostics, but the pressure projection still uses constant `rho`.

### `enable_variable_nu`

This controls whether Cantera viscosity may affect the flow viscosity. If `.true.`, it requires:

```fortran
enable_cantera = .true.
```

For current non-reacting scalar validation, the most useful setting is often:

```fortran
&fluid_input
  enable_cantera = .false.
  enable_variable_nu = .false.
/

&species_input
  enable_cantera = .true.   ! Cantera D_k only
/
```

This keeps the flow Reynolds number fixed while allowing species diffusivity to vary with `T`, `Y`, and `p0`.

### `background_temp`

Fallback temperature for Cantera transport when energy is disabled.

Current behavior:

```text
if energy is enabled:  Cantera transport uses energy%T
if energy is disabled: Cantera transport uses background_temp
```

### `background_press`

Uniform thermodynamic pressure used by Cantera.

Current behavior:

```text
p0 = background_press
```

Cantera uses this pressure for:

- `h(T,Y,p0)`
- `T(h,Y,p0)`
- `cp(T,Y,p0)`
- `thermal_conductivity(T,Y,p0)`
- `rho_thermo(T,Y,p0)`
- `D_k(T,Y,p0)`

Do not confuse this with the hydrodynamic/projection pressure field.

### `transport_update_interval`

Controls how often Cantera transport properties are refreshed.

It controls:

- Cantera `mu`, if variable viscosity is enabled.
- Cantera species diffusivity `D_k`, if species Cantera is enabled.

It does not control:

- Energy-side Cantera `h <-> T` recovery.
- `cp` refresh.
- Energy-side `thermal_conductivity` refresh.
- Projection density.
- Pressure projection.

Recommended values:

```text
1       debugging and validation
5-20    faster non-reacting development runs
```

## 5. Energy configuration

The `&energy_input` block controls enthalpy/temperature storage, energy transport, constant-property fallback values, and Cantera thermodynamics.

```fortran
&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.

  thermo_update_interval = 1
  thermo_default_species = "N2"

  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/
```

### `enable_energy`

Master switch for energy transport.

When enabled, the solver stores and advances:

- `T`: temperature `[K]`
- `h`: transported sensible enthalpy `[J/kg]`
- `qrad`: volumetric energy source `[W/m^3]`
- `cp`: heat capacity `[J/kg/K]`
- `lambda`: thermal conductivity `[W/m/K]`
- `rho_thermo`: Cantera diagnostic density `[kg/m^3]`

Current energy discretization:

```text
energy advection: upwind
energy diffusion: central/two-point gradient of T
source: qrad / rho
integration: explicit
```

Important: even though `h` is transported, conduction is driven by `grad(T)`.

### `enable_cantera_thermo`

When `.false.`, the solver uses the constant-cp fallback:

```text
h = h_ref + cp * (T - T_ref)
T = T_ref + (h - h_ref) / cp
lambda = energy_lambda
cp = energy_cp
```

When `.true.`, the solver uses Cantera thermodynamics:

```text
h_sensible = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
T = Cantera_HP_inverse(h_sensible + h_abs(T_ref,Y,p0), Y, p0)
cp = cp(T,Y,p0)
lambda = thermal_conductivity(T,Y,p0)
rho_thermo = rho(T,Y,p0)
```

Why sensible enthalpy matters:

Cantera absolute enthalpy includes formation/reference contributions. In non-reacting mixing, raw absolute enthalpy can create artificial heat release when composition changes. Using sensible enthalpy relative to the same composition at `T_ref` removes that artifact.


### Cantera thermo sync optimization

During the energy step, the solver may use a combined thermo-sync call:

```text
(T, cp, lambda, rho_thermo) = sync(h, Y, p0)
```

This is logically equivalent to recovering `T` from `h,Y,p0` and then refreshing `cp`, `lambda`, and `rho_thermo`, but avoids a redundant second Cantera pass.

For species-enabled runs, a pre-flux sync is used after species transport because `Y` may have changed. For species-disabled runs, the pre-flux sync can be skipped after initialization because the previous post-flux sync remains valid.

### Option A: species-energy coupling convention

This is now the accepted convention.

When species transport changes the composition from `Y_old` to `Y_new`, preserve transported enthalpy and recover temperature:

```text
h_after_species = h_before_species
T_after_species = T(h_after_species, Y_new, p0)
```

Do not recompute enthalpy from the old temperature and new composition before energy transport:

```text
Do not use: h = h(T_old, Y_new, p0)
```

That would preserve temperature instead of the transported energy state and can numerically add or remove sensible enthalpy during passive species mixing.

### `thermo_update_interval`

Reserved for future optimization. Current supported value:

```fortran
thermo_update_interval = 1
```

Cantera energy thermodynamics must remain logically synchronized every energy step. The implementation may optimize internally with a combined thermo-sync call and conservative cache, but do not set this above 1 until a stale-thermo strategy is explicitly implemented and validated.

### `thermo_default_species`

Default species for Cantera thermodynamics when species transport is off.

Example:

```fortran
thermo_default_species = "N2"
```

Used when:

```text
enable_cantera_thermo = .true.
enable_species = .false.
```

Ignored when:

```text
enable_species = .true.
```

In that case, Cantera thermodynamics use transported `species%Y`.

Useful examples:

```fortran
thermo_default_species = "N2"
thermo_default_species = "CO2"
thermo_default_species = "H2O"
```

Use `"H2O"` with the letter `O`, not `"H20"` with zero.

### `initial_T`

Initial gas temperature in Kelvin. Must be positive.

### `energy_reference_T`

Reference temperature used for sensible enthalpy. Must be positive.

Recommended for simple validation:

```fortran
initial_T = 300.0
energy_reference_T = 300.0
```

This makes the initial sensible enthalpy close to zero for the initial composition.

### `energy_reference_h`

Reference enthalpy for constant-cp mode. Usually:

```fortran
energy_reference_h = 0.0
```

### `energy_cp`

Constant heat capacity for non-Cantera energy mode. Must be positive.

### `energy_lambda`

Constant thermal conductivity for non-Cantera energy mode. Must be non-negative.

For visible diffusion debugging, you may temporarily use a larger value than air's physical value.

## 6. Species configuration

The `&species_input` block controls passive species transport, fallback diffusivity, Cantera species diffusivity, and future reactions.

```fortran
&species_input
  enable_species = .true.
  enable_reactions = .false.
  enable_cantera = .true.

  nspecies = 3
  species_name = "O2", "N2", "CO2"

  initial_Y = 0.0, 1.0, 0.0
  species_diffusivity = 2.0e-5, 2.0e-5, 2.0e-5
/
```

### `enable_species`

Master switch for species transport.

### `enable_reactions`

Future reaction chemistry switch. Current validation setting:

```fortran
enable_reactions = .false.
```

Do not turn reactions on until passive species, passive energy, and Cantera thermo have been validated.

### `enable_cantera` in `species_input`

This controls Cantera species diffusivity `D_k`.

When `.false.`:

```text
D_k = species_diffusivity(k)
```

When `.true.`:

```text
D_k = Cantera mixture-averaged diffusivity D_k(T,Y,p0)
```

The temperature source is:

```text
energy enabled:  energy%T
energy disabled: background_temp
```

### `nspecies`

Number of transported species. Must be non-negative and within the compiled `max_species` limit.

### `species_name`

Transported species names. These must match species in the Cantera mechanism whenever Cantera species transport or Cantera thermodynamics with transported species is enabled.

Example:

```fortran
species_name = "O2", "N2", "CO2"
```

### `initial_Y`

Initial domain mass fractions.

Examples:

```fortran
initial_Y = 0.0, 1.0, 0.0   ! pure N2 for O2,N2,CO2 list
initial_Y = 0.0, 0.0, 1.0   ! pure CO2 for O2,N2,CO2 list
```

### `species_diffusivity`

Fallback constant diffusivity values `[m^2/s]`. Used when `species_input enable_cantera = .false.`.

## 7. Boundary configuration

Boundary configuration links the physical boundary names in `patches.dat` to flow, pressure, species, and temperature conditions.

Example:

```fortran
&boundary_input
  n_patches = 5
  patch_name = "inlet", "outlet", "wall", "zmin", "zmax"

  patch_type = "dirichlet", "neumann", "wall", "symmetric", "symmetric"

  patch_velocity_type = "fixed_value", "zero_gradient", "no_slip", "symmetric", "symmetric"
  patch_pressure_type = "zero_gradient", "zero_gradient", "zero_gradient", "zero_gradient", "zero_gradient"
  patch_temperature_type = "fixed_value", "zero_gradient", "zero_gradient", "zero_gradient", "zero_gradient"
  patch_species_type = "fixed_value", "zero_gradient", "zero_gradient", "zero_gradient", "zero_gradient"

  patch_u = 1.0, 0.0, 0.0, 0.0, 0.0
  patch_v = 0.0, 0.0, 0.0, 0.0, 0.0
  patch_w = 0.0, 0.0, 0.0, 0.0, 0.0

  patch_p = 0.0, 0.0, 0.0, 0.0, 0.0
  patch_dpdn = 0.0, 0.0, 0.0, 0.0, 0.0

  patch_T = 400.0, 300.0, 300.0, 300.0, 300.0

  patch_Y(1,1) = 0.0
  patch_Y(2,1) = 0.0
  patch_Y(3,1) = 1.0
/
```

### Boundary type aliases

The parser recognizes these aliases:

| Physical meaning | Accepted names |
|---|---|
| Wall | `"wall"`, `"no_slip"`, `"moving_wall"` |
| Symmetry/slip | `"symmetry"`, `"symmetric"`, `"slip"` |
| Periodic | `"periodic"` |
| Fixed value | `"dirichlet"`, `"fixed_value"` |
| Zero gradient | `"neumann"`, `"zero_gradient"` |

### Velocity boundaries

- Wall/no-slip uses the configured wall velocity.
- Symmetry removes the normal velocity component.
- Zero-gradient and periodic use extrapolated/effective neighbor behavior.

### Pressure boundaries

Use fixed pressure only where needed to anchor the pressure field. Otherwise use zero-gradient for typical outlets/walls, depending on the case design.

### Temperature boundaries

Use:

```fortran
patch_temperature_type = "fixed_value"
patch_T = ...
```

for fixed-temperature inlets or hot/cold walls.

If a patch has:

```fortran
patch_temperature_type = "zero_gradient"
```

then `patch_T` is not used as a fixed value on that patch.

### Species boundaries

Use:

```fortran
patch_species_type = "fixed_value"
patch_Y(k,patch) = ...
```

for species inlets.

Mass fractions should sum to about 1 on fixed-value species patches.

### Fixed-temperature inlet with species composition

For Cantera thermo plus transported species, a fixed-temperature species inlet should be interpreted as:

```text
T_b = patch_T(patch)
Y_b = patch_Y(:,patch)
h_b = h(T_b,Y_b,p0)
```

This is important for hot/cold fuel/oxidizer inlets. It avoids evaluating boundary enthalpy with interior composition or a default bath gas.

## 8. Output configuration

```fortran
&output_input
  output_dir = "cases/rectangle_2D/output"
  write_vtu = .true.
  write_diagnostics = .true.
/
```

### `output_dir`

Directory for generated output files.

Before rerunning after changing output arrays, remove old output:

```bash
rm -rf cases/rectangle_2D/output
```

### `write_vtu`

Controls VTU/PVTU/PVD visualization output.

### `write_diagnostics`

Controls CSV diagnostics. Current diagnostic files may include:

- `diagnostics.csv`
- `energy_diagnostics.csv`

## 9. Profiling configuration

```fortran
&profiling_input
  enable_profiling = .true.
  nested_profiling = .true.
/
```

### `enable_profiling`

Enable MPI-aware profiling.

### `nested_profiling`

Show a nested call tree in the terminal profiling report.


Current profiling region names include:

```text
Transport_Update
Projection_Step
Species_Transport
Energy_Transport
Diagnostics_Write_Flow
Diagnostics_Write_Energy
Output_Write_VTU
```

For Cantera thermo energy runs, nested energy children may include:

```text
Energy_Cantera_PreSync
Energy_PreFlux_Exchange
Energy_Flux_Update
Energy_Cantera_PostSync
Energy_Final_Exchange
```

The profiler reports inclusive wall time. Flat rows are not additive when nested timers are enabled.

## 10. Recommended configuration modes

### Mode 1: Baseline hydrodynamics

Use this to verify that energy/species additions did not change old flow behavior.

```fortran
&fluid_input
  rho = 1.0
  nu = 1.0e-2
  enable_cantera = .false.
  enable_variable_density = .false.
  enable_variable_nu = .false.
  cantera_mech_file = "gri30.yaml"
  background_temp = 300.0
  background_press = 101325.0
  transport_update_interval = 10
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
  thermo_update_interval = 1
  thermo_default_species = "N2"
  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/
```

Expected result:

- Old flow diagnostics remain unchanged within numerical tolerance.
- No energy fields are required.

### Mode 2: Constant-property energy

Use this to validate energy numerics without Cantera complexity.

```fortran
&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .false.
  thermo_update_interval = 1
  thermo_default_species = "N2"
  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/
```

Expected result:

- Uniform temperature remains uniform.
- Pure diffusion is smooth and stable under explicit timestep limits.
- `qrad` remains zero unless a source is prescribed later.

### Mode 3: Cantera thermodynamics without transported species

Use this for Stage 2A thermo smoke tests.

```fortran
&species_input
  enable_species = .false.
  enable_reactions = .false.
  enable_cantera = .false.
  nspecies = 0
/

&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.
  thermo_update_interval = 1
  thermo_default_species = "N2"
  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/
```

Change `thermo_default_species` from `"N2"` to `"CO2"` to verify that `cp`, `thermal_conductivity`, `rho_thermo`, and thermal diffusivity respond to the selected gas.

### Mode 4: Cantera thermodynamics with transported species

Use this for Stage 2B non-reacting species/energy coupling.

```fortran
&species_input
  enable_species = .true.
  enable_reactions = .false.
  enable_cantera = .false.
  nspecies = 3
  species_name = "O2", "N2", "CO2"
  initial_Y = 0.0, 1.0, 0.0
  species_diffusivity = 2.0e-5, 2.0e-5, 2.0e-5
/

&energy_input
  enable_energy = .true.
  enable_cantera_thermo = .true.
  thermo_update_interval = 1
  thermo_default_species = "N2"
  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/
```

In this mode:

```text
thermo_default_species is ignored
Cantera thermo uses transported species%Y
h is preserved through species changes
T is recovered from h,Y,p0
```

### Mode 5: Cantera species diffusivity with fixed-Re flow

Use this to test variable `D_k` while keeping the flow Reynolds number fixed.

```fortran
&fluid_input
  rho = 1.0
  nu = 1.0e-2
  enable_cantera = .false.
  enable_variable_density = .false.
  enable_variable_nu = .false.
  cantera_mech_file = "gri30.yaml"
  background_temp = 300.0
  background_press = 101325.0
  transport_update_interval = 10
/

&species_input
  enable_species = .true.
  enable_reactions = .false.
  enable_cantera = .true.
  nspecies = 3
  species_name = "O2", "N2", "CO2"
  initial_Y = 0.0, 1.0, 0.0
  species_diffusivity = 2.0e-5, 2.0e-5, 2.0e-5
/
```

Expected result:

```text
D_k = Cantera D_k(T,Y,p0)
nu = constant
rho = constant
Re = fixed
```

### Mode 6: Non-reacting counterflow development target

Use this only after simpler tests pass.

Recommended rules:

- `enable_species = .true.`
- `enable_reactions = .false.`
- `enable_energy = .true.`
- `enable_cantera_thermo = .true.`
- `enable_variable_density = .false.`
- `enable_variable_nu = .false.` initially
- Fixed-temperature/fixed-species inlets use matching `patch_T` and `patch_Y`
- Use Option A enthalpy/species coupling

Expected result:

- Non-reacting mixing layer trends are plausible.
- Temperature remains physical and bounded.
- `sum_Y` remains close to 1.
- Exact agreement with fully variable-density counterflow is not expected yet.

## 11. ParaView and VTU fields

Depending on enabled physics, ParaView may show the following arrays.

### Flow fields

- `velocity`
- `pressure`

`pressure` is hydrodynamic/projection pressure, not thermodynamic pressure.

### Energy fields

- `temperature`
- `enthalpy`
- `qrad`
- `cp`
- `thermal_conductivity`
- `rho_thermo`
- `thermal_diffusivity`
- `thermo_pressure`

Meanings:

| Field | Meaning |
|---|---|
| `temperature` | Active gas temperature `[K]`. |
| `enthalpy` | Transported sensible enthalpy `[J/kg]`. |
| `qrad` | Volumetric radiation/source term `[W/m^3]`. Positive adds energy to gas. |
| `cp` | Constant or Cantera heat capacity `[J/kg/K]`. |
| `thermal_conductivity` | Constant or Cantera conductivity `[W/m/K]`. |
| `rho_thermo` | Cantera diagnostic density `[kg/m^3]`. Not used in projection. |
| `thermal_diffusivity` | `thermal_conductivity / (rho_flow * cp)`. |
| `thermo_pressure` | Uniform `background_press`. |

### Species fields

Examples:

- `Y_O2`
- `Y_N2`
- `Y_CO2`
- `sum_Y`
- `D_O2`
- `D_N2`
- `D_CO2`

`D_<species>` is the diffusivity used by species transport.

If arrays are missing in ParaView after adding new output fields, remove old output and rerun:

```bash
rm -rf cases/rectangle_2D/output
make rectangle_2D-release NP=8
```

Then inspect the output metadata:

```bash
grep -R "thermal_diffusivity\|thermo_pressure\|D_" cases/rectangle_2D/output
```

## 12. Runtime terminal reporting

The solver should print a mode summary similar to:

```text
Flow density mode: constant rho = ... kg/m^3
Flow viscosity mode: constant nu = ... m^2/s
Cantera rho_thermo: diagnostic only, not used by projection
Cantera transport update interval: N step(s) [mu/D_k only]
Cantera transport temperature source: energy%T
Cantera thermo update interval: every energy step
Cantera thermo sync: combined T/cp/lambda/rho_thermo from h,Y,p0
Cantera thermodynamic pressure p0: ... Pa
```

Use this block to verify that the active runtime mode matches the intended validation mode.

## 13. Recommended validation ladder

### Step 1: Energy disabled regression

Run an existing flow-only or flow/species case with:

```fortran
enable_energy = .false.
```

Expected:

- Old flow/species behavior is unchanged.
- Build and output still work.

### Step 2: Constant-property energy initialization

Use:

```fortran
enable_energy = .true.
enable_cantera_thermo = .false.
```

Expected:

- `temperature` is initialized to `initial_T`.
- `enthalpy` follows constant-cp formula.
- `qrad = 0`.

### Step 3: Constant-property diffusion

Use zero velocity and fixed hot/cold temperature boundaries.

Expected:

- Temperature evolves smoothly.
- Steady profile is reasonable for the geometry.
- Explicit timestep is stable.

### Step 4: Cantera h/T roundtrip without species

Use:

```fortran
enable_cantera_thermo = .true.
enable_species = .false.
```

Expected:

- `h(T)` followed by `T(h)` recovers the original temperature.
- Changing `thermo_default_species` changes thermodynamic properties.

### Step 5: Cantera thermo with species, no reactions

Use:

```fortran
enable_species = .true.
enable_reactions = .false.
enable_cantera_thermo = .true.
```

Expected:

- Thermodynamic properties follow transported composition.
- Enthalpy is preserved through species updates using Option A.
- Temperature remains bounded during passive mixing.

### Step 6: Cantera species diffusivity with fixed-Re flow

Use species Cantera only:

```fortran
&fluid_input
  enable_cantera = .false.
  enable_variable_nu = .false.
/

&species_input
  enable_cantera = .true.
/
```

Expected:

- `D_k` varies with `T`, `Y`, and `p0`.
- Flow viscosity remains fixed.

### Step 7: Non-reacting counterflow

Use opposed inlets, species boundary compositions, fixed inlet temperatures, and reactions off.

Compare:

- Centerline temperature.
- Centerline species profiles.
- Mixture fraction trend.
- Stagnation-plane location.
- Scalar dissipation trend.

Expected:

- Qualitative non-reacting mixing-layer behavior.
- No expectation of exact variable-density Cantera 1D agreement until variable-density low-Mach coupling is implemented.

### Step 8: Manufactured `qrad` tests, future

Before external radiation physics, test prescribed `qrad` patterns:

- Uniform heating.
- Uniform cooling.
- Gaussian heating/cooling.
- Localized source.

Expected:

- Positive `qrad` heats the gas.
- Negative `qrad` cools the gas.
- Domain-integrated energy change matches integrated `qrad`.

## 14. Common mistakes

### Mistake 1: Enabling variable density too early

Do not use:

```fortran
enable_variable_density = .true.
```

Current solver does not implement variable-density projection. Use `rho_thermo` only as a diagnostic field.

### Mistake 2: Turning on all Cantera switches at once

Avoid starting with:

```fortran
fluid_input enable_cantera = .true.
species_input enable_cantera = .true.
enable_cantera_thermo = .true.
enable_variable_nu = .true.
```

Instead, validate one feature at a time.

### Mistake 3: Confusing the two `enable_cantera` switches

There are two separate blocks:

```fortran
&fluid_input
  enable_cantera = ...   ! flow viscosity path
/

&species_input
  enable_cantera = ...   ! species diffusivity path
/
```

Use the species switch for Cantera `D_k`. Use the fluid switch only when testing Cantera viscosity.

### Mistake 4: Using projection pressure as thermodynamic pressure

Current thermodynamic pressure is:

```fortran
background_press = 101325.0
```

The pressure field written to output is projection pressure, not Cantera thermodynamic pressure.

### Mistake 5: Fixed temperature with missing temperature type

Setting `patch_T` alone is not enough. You also need:

```fortran
patch_temperature_type = "fixed_value"
```

for that patch.

### Mistake 6: Species inlet mass fractions not summing to 1

For fixed species inlets, check:

```text
sum_k patch_Y(k,patch) ~= 1
```

This is especially important for Cantera thermo and species diffusivity.

### Mistake 7: Using interior composition for inlet enthalpy

For fixed-temperature species inlets, the intended rule is:

```text
h_b = h(T_b,Y_b,p0)
```

not:

```text
h_b = h(T_b,Y_interior,p0)
```

### Mistake 8: Recomputing enthalpy after species transport

The accepted Option A convention is:

```text
preserve h
recover T from h,Y_new,p0
```

Do not preserve old temperature by recomputing `h(T_old,Y_new,p0)`.

### Mistake 9: Expecting exact counterflow agreement too early

The current flow solver remains constant-density. Exact agreement with variable-density 1D Cantera counterflow is not expected yet.

### Mistake 10: Reusing old ParaView output after changing arrays

Remove output when array structure changes:

```bash
rm -rf cases/<case>/output
```

Then rerun the case.

## 15. Suggested documentation set in the repo

Recommended files:

```text
docs/case_nml_reference.md
docs/input_configuration_guide.md
docs/energy_thermo_conventions.md
docs/enthalpy_species_coupling_convention.md
docs/validation_ladder.md
```

This file can serve as `docs/input_configuration_guide.md`. The compact companion file can serve as `docs/case_nml_reference.md`.
