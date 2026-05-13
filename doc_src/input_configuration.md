title: case.nml Reference

# LowMachReact-Hex `case.nml` Reference

This file is the compact reference for writing and reviewing LowMachReact-Hex simulation input files. It lists the expected namelist blocks, available options, valid values, and the current meaning of each switch.

Use this file when you want to quickly check the structure of `case.nml`. Use `input_configuration_guide.md` for longer explanations, validation workflows, and examples.

## Current solver stage

The current development state is:

- Flow projection is constant-density.
- `rho` is the flow/projection density.
- `rho_thermo` is a Cantera diagnostic field only.
- Energy transport is available.
- The transported energy variable is sensible enthalpy `h`.
- Temperature `T` is recovered from `h`, composition `Y`, and thermodynamic pressure `p0`.
- Cantera thermodynamics can provide `h <-> T`, `cp`, `thermal_conductivity`, and `rho_thermo`.
- Cantera species diffusivity `D_k` can be enabled independently from Cantera flow viscosity.
- Real radiation physics is not coupled yet. `qrad` exists as the volumetric source interface.
- Option A convention is used for species-energy coupling: preserve transported `h`, then recover `T = T(h,Y,p0)` after composition changes.
- Fixed-temperature boundaries with species enabled should use boundary composition: `h_boundary = h(T_boundary,Y_boundary,p0)`.

## Canonical file skeleton

```fortran
&mesh_input
  mesh_dir = "cases/rectangle_2D/mesh_native"
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
  enable_variable_nu = .false.
  cantera_mech_file = "gri30.yaml"
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
  thermo_update_interval = 1
  thermo_default_species = "N2"
  initial_T = 300.0
  energy_reference_T = 300.0
  energy_reference_h = 0.0
  energy_cp = 1005.0
  energy_lambda = 0.026
/

&output_input
  output_dir = "output"
  write_vtu = .true.
  write_diagnostics = .true.
/

&profiling_input
  enable_profiling = .false.
  nested_profiling = .false.
/
```

Notes:

- Fortran logical values are `.true.` and `.false.`.
- Quote strings.
- Namelist blocks are read by name; keeping the order above is recommended for readability.
- Unknown variables in strict blocks should be treated as input errors.

## `&mesh_input`

| Option | Type | Default/example | Description |
|---|---:|---|---|
| `mesh_dir` | string | `"cases/.../mesh_native"` | Directory containing native mesh files. |

Required files in `mesh_dir`:

- `points.dat`
- `cells.dat`
- `faces.dat`
- `patches.dat`
- `periodic.dat` only when periodic boundaries are used

## `&time_input`

| Option | Type | Default/example | Description |
|---|---:|---|---|
| `nsteps` | integer | `1000` | Total number of timesteps. Must be non-negative. |
| `dt` | real | `1.0e-4` | Fixed timestep, or initial timestep when dynamic CFL control is enabled. Must be positive. |
| `output_interval` | integer | `100` | Output and diagnostics cadence in steps. Must be positive. |
| `use_dynamic_dt` | logical | `.false.` | Enables CFL-based timestep adjustment. |
| `max_cfl` | real | `0.5` | Target CFL when dynamic timestep is enabled. |

## `&fluid_input`

| Option | Type | Default/example | Description |
|---|---:|---|---|
| `rho` | real | `1.0` | Constant flow/projection density `[kg/m^3]`. Used by the current projection solver. |
| `nu` | real | `1.0e-2` | Constant kinematic viscosity `[m^2/s]`. Used when variable Cantera viscosity is off. |
| `enable_cantera` | logical | `.false.` | Enables Cantera fluid transport path for viscosity when paired with `enable_variable_nu`. This is separate from `species_input enable_cantera`. |
| `enable_variable_density` | logical | `.false.` | Reserved future flag. Current supported value is `.false.` only. `.true.` should stop because variable-density projection is not implemented. |
| `enable_variable_nu` | logical | `.false.` | If `.true.`, Cantera viscosity may affect flow viscosity. Requires `fluid_input enable_cantera = .true.`. |
| `cantera_mech_file` | string | `"gri30.yaml"` | Cantera mechanism path. Required for Cantera transport or thermo modes. |
| `background_temp` | real | `300.0` | Fallback transport temperature when energy is disabled. |
| `background_press` | real | `101325.0` | Uniform thermodynamic pressure `p0` used by Cantera. Not projection pressure. |
| `transport_update_interval` | integer | `10` | Refresh interval for Cantera transport properties `mu` and/or `D_k`. Must be positive. |

Recommended current validation settings:

```fortran
rho = 1.0
nu = 1.0e-2
enable_variable_density = .false.
enable_variable_nu = .false.
```

Use `species_input enable_cantera = .true.` with `fluid_input enable_cantera = .false.` when you want variable Cantera species diffusivity but fixed Reynolds number.

## `&solver_input`

| Option | Type | Default/example | Description |
|---|---:|---|---|
| `pressure_max_iter` | integer | `2000` | Maximum pressure Poisson iterations. Must be positive. |
| `pressure_tol` | real | `1.0e-8` | Pressure solver tolerance. Must be positive. |
| `body_force_x` | real | `0.0` | Body force x component. |
| `body_force_y` | real | `0.0` | Body force y component. |
| `body_force_z` | real | `0.0` | Body force z component. |
| `convection_scheme` | string | `"central"` | Flow/momentum convection scheme. Common values: `"central"`, `"upwind"`. |

Current note: `convection_scheme` controls the flow/momentum path. Energy advection currently uses upwind for robustness. Higher-order bounded scalar schemes should be a future numerics stage.

## `&boundary_input`

| Option | Type | Description |
|---|---:|---|
| `n_patches` | integer | Number of mesh patches. Must match `patches.dat`. |
| `patch_name(max_patches)` | string array | Patch names matching the mesh. |
| `patch_type(max_patches)` | string array | General patch type. |
| `patch_velocity_type(max_patches)` | string array | Velocity-specific boundary type. |
| `patch_pressure_type(max_patches)` | string array | Pressure-specific boundary type. |
| `patch_temperature_type(max_patches)` | string array | Temperature/enthalpy boundary type. |
| `patch_species_type(max_patches)` | string array | Species boundary type. |
| `patch_u(max_patches)` | real array | Boundary velocity x value. |
| `patch_v(max_patches)` | real array | Boundary velocity y value. |
| `patch_w(max_patches)` | real array | Boundary velocity z value. |
| `patch_p(max_patches)` | real array | Boundary pressure value for fixed pressure boundaries. |
| `patch_dpdn(max_patches)` | real array | Boundary pressure gradient value for gradient boundaries. |
| `patch_T(max_patches)` | real array | Boundary temperature `[K]`. Used only when temperature type is fixed value. |
| `patch_Y(species_id,patch_id)` | real array | Species boundary mass fractions. Used only for fixed species boundaries. |

Recognized boundary aliases:

| Meaning | Accepted names |
|---|---|
| Wall | `"wall"`, `"no_slip"`, `"moving_wall"` |
| Symmetry/slip | `"symmetry"`, `"symmetric"`, `"slip"` |
| Periodic | `"periodic"` |
| Dirichlet/fixed value | `"dirichlet"`, `"fixed_value"` |
| Neumann/zero gradient | `"neumann"`, `"zero_gradient"` |

Temperature boundary rule:

- `patch_temperature_type = "fixed_value"` uses `patch_T`.
- Other temperature types behave as zero-gradient for temperature at this stage.
- With Cantera thermo and species enabled, fixed-temperature inlet enthalpy should be evaluated with boundary composition: `h(T_b,Y_b,p0)`.

Species boundary rule:

- `patch_species_type = "fixed_value"` uses `patch_Y(:,patch)`.
- For fixed species boundaries, mass fractions should sum to about 1.

## `&species_input`

| Option | Type | Default/example | Description |
|---|---:|---|---|
| `enable_species` | logical | `.false.` | Master switch for species transport. |
| `enable_reactions` | logical | `.false.` | Reserved/future reaction switch. Keep `.false.` for current non-reacting validation. |
| `enable_cantera` | logical | `.false.` | Enables Cantera species diffusivity `D_k`. This is separate from `fluid_input enable_cantera`. |
| `nspecies` | integer | `3` | Number of transported species. Must be within compiled limit. |
| `species_name(max_species)` | string array | `"O2", "N2", "CO2"` | Names of transported species. Must exist in Cantera mechanism when Cantera is used. |
| `initial_Y(max_species)` | real array | `0.0, 1.0, 0.0` | Initial mass fractions. |
| `species_diffusivity(max_species)` | real array | `2.0e-5` | Constant fallback diffusivities `[m^2/s]`. Used when species Cantera is off. |

Current recommended non-reacting validation:

```fortran
enable_species = .true.
enable_reactions = .false.
enable_cantera = .true.   ! variable D_k only, if desired
```

## `&energy_input`

| Option | Type | Default/example | Description |
|---|---:|---|---|
| `enable_energy` | logical | `.false.` | Master switch for enthalpy/temperature energy transport. |
| `enable_cantera_thermo` | logical | `.false.` | Enables Cantera `h <-> T`, `cp`, `lambda`, and diagnostic `rho_thermo`. |
| `thermo_update_interval` | integer | `1` | Reserved optimization knob. Current supported value is `1`; Cantera thermo remains logically synchronized every energy step. |
| `thermo_default_species` | string | `"N2"` | Default species used for Cantera thermo when species transport is off. Ignored when transported species are enabled. |
| `initial_T` | real | `300.0` | Initial temperature `[K]`. Must be positive. |
| `energy_reference_T` | real | `300.0` | Reference temperature for sensible enthalpy `[K]`. Must be positive. |
| `energy_reference_h` | real | `0.0` | Reference enthalpy for constant-cp fallback `[J/kg]`. |
| `energy_cp` | real | `1005.0` | Constant heat capacity `[J/kg/K]` when Cantera thermo is off. Must be positive. |
| `energy_lambda` | real | `0.026` | Constant thermal conductivity `[W/m/K]` when Cantera thermo is off. Must be non-negative. |

Current energy convention:

```text
transported state: h_sensible [J/kg]
dependent state:   T [K]
recovery:          T = T(h,Y,p0)
p0:                background_press
```

Cantera sensible enthalpy convention:

```text
h_sensible(T,Y,p0) = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
```


Combined Cantera thermo sync used by the energy step:

```text
(T, cp, lambda, rho_thermo) = sync(h_sensible,Y,p0)
```

The combined sync preserves transported `h_sensible` and may reuse cached dependent thermo state only when `h_sensible`, `Y`, and `p0` are unchanged within tight tolerances.

## `&output_input`

| Option | Type | Default/example | Description |
|---|---:|---|---|
| `output_dir` | string | `"output"` | Directory for output files. Must be non-empty. |
| `write_vtu` | logical | `.true.` | Enables visualization output. |
| `write_diagnostics` | logical | `.true.` | Enables CSV diagnostics. |

Possible outputs include:

- `diagnostics.csv`
- `energy_diagnostics.csv`
- VTU/PVTU/PVD files

## `&profiling_input`

| Option | Type | Default/example | Description |
|---|---:|---|---|
| `enable_profiling` | logical | `.false.` | Enable MPI-aware profiler. |
| `nested_profiling` | logical | `.false.` | Show nested profiler call tree. Current timer names include `Transport_Update`, `Projection_Step`, `Species_Transport`, `Energy_Transport`, and output/diagnostic split timers. |

## Minimal hydrodynamic case

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
  enable_variable_nu = .false.
  cantera_mech_file = "gri30.yaml"
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
/
```

Fill in `boundary_input` for real cases; the empty block is shown only to illustrate a minimal option layout.
