title: case.nml Reference


# LowMachReact-Hex `case.nml` Reference

This is the compact reference for writing and reviewing LowMachReact-Hex
simulation input files.

Use this file for quick option lookup.  Use `input_configuration_guide.md` for
longer explanations, validation workflows, examples, and common mistakes.

This reference reflects the current implementation:

```text
stable baseline:
  constant-density incompressible / low-Mach projection solver

active guarded path:
  non-reacting variable-density low-Mach solver using Cantera density,
  conservative rho*Y and rho*h transport, and low-Mach source diagnostics

future work:
  reactions, reaction heat release, physical radiation solver, restart, and
  broader EOS/pressure validation
```

---

## 1. Current mode summary

### 1.1 Constant-density mode

```fortran
enable_variable_density = .false.
density_eos = "constant"
```

Meaning:

```text
transport%rho = params%rho
projection target = div(u) = 0
fields%mass_flux = rho * fields%face_flux
energy%rho_thermo, if present, is diagnostic
```

### 1.2 Cantera-assisted constant-density mode

```fortran
enable_variable_density = .false.
density_eos = "constant"
enable_energy = .true.
enable_cantera_thermo = .true.
```

Meaning:

```text
flow/projection density = params%rho
Cantera provides h <-> T, cp, lambda, rho_thermo
rho_thermo is not active density
```

### 1.3 Guarded variable-density low-Mach mode

Supported active combination:

```fortran
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
thermo_update_interval = 1
enable_reactions = .false.
```

Meaning:

```text
transport%rho <- energy%rho_thermo
projection target = div(u) = S_projection
species branch = conservative rho*Y
energy branch = conservative rho*h
```

`density_eos="cantera"` means active density comes from the selected Cantera
phase.  The actual EOS is the selected phase's `thermo:` model.

### 1.4 Reserved or future modes

```text
density_eos = "ideal_gas"      parsed/reserved for future solver-side EOS
enable_reactions = .true.      future chemistry/heat-release path
full multicomponent FV fluxes  future species-operator extension
physical radiation solver      future qrad producer
fully compressible flow        outside current architecture
```

---

## 2. Canonical file skeleton

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
  output_dir = "output"
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

Notes:

```text
- Strings must be quoted.
- Logical values are .true. and .false.
- Namelist blocks are read by name; the order above is recommended.
- Unknown variables should be treated as input errors.
- Fill in boundary_input for real cases.  n_patches=0 is only a skeleton placeholder.
```

---

## 3. `&mesh_input`

| Option | Type | Required | Description |
|---|---:|---:|---|
| `mesh_dir` | string | yes | Directory containing native mesh files. |

Required files:

```text
points.dat
cells.dat
faces.dat
patches.dat
periodic.dat  only when periodic boundaries are used
```

Example:

```fortran
&mesh_input
  mesh_dir = "cases/channel_flow/mesh_native"
/
```

---

## 4. `&time_input`

| Option | Type | Required | Description |
|---|---:|---:|---|
| `nsteps` | integer | yes | Final timestep count. Must be non-negative. |
| `dt` | real | yes | Fixed timestep, or initial timestep if dynamic CFL mode is enabled. Must be positive. |
| `output_interval` | integer | yes | Output and diagnostics cadence in timesteps. Must be positive. |
| `use_dynamic_dt` | logical | no | Enable CFL-based timestep adjustment. |
| `max_cfl` | real | if dynamic dt | Target CFL when `use_dynamic_dt=.true.`. |

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

---

## 5. `&fluid_input`

| Option | Type | Description |
|---|---:|---|
| `rho` | real | Constant/reference flow density `[kg/m^3]`. Active density in constant-density mode; startup/reference value in variable-density mode. |
| `nu` | real | Constant kinematic viscosity `[m^2/s]` when variable viscosity is off. |
| `enable_cantera` | logical | Enables Cantera fluid transport path, mainly viscosity. Separate from `species_input enable_cantera`. |
| `enable_variable_density` | logical | Enables guarded variable-density low-Mach path when used with supported settings. |
| `density_eos` | string | Active-density selector: `"constant"`, `"cantera"`, or reserved `"ideal_gas"`. |
| `enable_variable_nu` | logical | If `.true.`, update `nu` from Cantera viscosity and active density. Requires fluid Cantera path. |
| `cantera_mech_file` | string | Cantera mechanism path or Cantera-resolvable mechanism name. |
| `cantera_phase_name` | string | Optional Cantera phase selector. Blank uses Cantera default/first phase. |
| `background_temp` | real | Fallback transport temperature `[K]` when energy is disabled. |
| `background_press` | real | Uniform thermodynamic pressure `p0` `[Pa]` for Cantera. Not projection pressure. |
| `transport_update_interval` | integer | Refresh interval for Cantera transport properties such as `mu` and `D_k`. Must be positive. |

### Valid `density_eos` values

| Value | Current meaning |
|---|---|
| `"constant"` | Use configured `rho` as active flow density. |
| `"cantera"` | Use selected Cantera phase density as active density only in guarded variable-density mode. |
| `"ideal_gas"` | Reserved for future solver-side ideal-gas EOS. |

### Cantera phase rule

```fortran
cantera_mech_file = "mechanisms/pr_high.yaml"
cantera_phase_name = "gas"
```

The selected YAML phase defines the actual EOS and transport model, for example:

```text
thermo: ideal-gas
thermo: Peng-Robinson
thermo: Redlich-Kwong

transport: mixture-averaged
transport: multicomponent
transport: high-pressure
transport: high-pressure-Chung
```

### Common fluid configurations

Fixed-Re constant-density validation:

```fortran
enable_cantera = .false.
enable_variable_density = .false.
density_eos = "constant"
enable_variable_nu = .false.
```

Variable-property non-reacting low-Mach:

```fortran
enable_cantera = .true.
enable_variable_density = .true.
density_eos = "cantera"
enable_variable_nu = .true.
```

---

## 6. `&solver_input`

| Option | Type | Description |
|---|---:|---|
| `pressure_max_iter` | integer | Maximum pressure Poisson iterations. Must be positive. |
| `pressure_tol` | real | Pressure solver tolerance. Must be positive. |
| `body_force_x` | real | Body force x component. |
| `body_force_y` | real | Body force y component. |
| `body_force_z` | real | Body force z component. |
| `convection_scheme` | string | Convection selector. Common values: `"central"`, `"upwind"`. |

Notes:

```text
- Use "upwind" for first stability tests.
- Use "central" for cleaner laminar validation after stability is established.
- Energy/species transport may use robust upwind paths internally even when flow
  convection is configured separately.
```

---

## 7. `&boundary_input`

| Option | Type | Description |
|---|---:|---|
| `n_patches` | integer | Number of configured mesh patches. |
| `patch_name(max_patches)` | string array | Patch names matching `patches.dat`. |
| `patch_type(max_patches)` | string array | Legacy/fallback patch type. Keep for compatibility. |
| `patch_velocity_type(max_patches)` | string array | Velocity-specific boundary type. |
| `patch_pressure_type(max_patches)` | string array | Pressure-specific boundary type. |
| `patch_temperature_type(max_patches)` | string array | Temperature/enthalpy boundary type. |
| `patch_species_type(max_patches)` | string array | Species boundary type. |
| `patch_u(max_patches)` | real array | Boundary velocity x component. |
| `patch_v(max_patches)` | real array | Boundary velocity y component. |
| `patch_w(max_patches)` | real array | Boundary velocity z component. |
| `patch_p(max_patches)` | real array | Fixed pressure value for pressure Dirichlet patches. |
| `patch_dpdn(max_patches)` | real array | Pressure normal gradient for Neumann patches. |
| `patch_T(max_patches)` | real array | Boundary temperature `[K]` for fixed-temperature patches. |
| `patch_Y(species_id,patch_id)` | real array | Boundary species mass fractions for fixed-species patches. |

### Boundary aliases

| Meaning | Accepted names |
|---|---|
| Wall | `"wall"`, `"no_slip"`, `"moving_wall"` |
| Symmetry/slip | `"symmetry"`, `"symmetric"`, `"slip"` |
| Periodic | `"periodic"` |
| Fixed value | `"dirichlet"`, `"fixed_value"` |
| Zero gradient | `"neumann"`, `"zero_gradient"` |

### Boundary rules

Velocity boundaries prescribe velocity and therefore volumetric flux:

```text
face_flux = u_b dot n * area
mass_flux = rho_face * face_flux
```

In variable-density mode, equal inlet velocity does not imply equal mass flux.

For fixed temperature:

```fortran
patch_temperature_type = "fixed_value"
patch_T = ...
```

For fixed composition:

```fortran
patch_species_type = "fixed_value"
patch_Y(k,p) = ...
```

For fixed-temperature plus fixed-composition inlets:

```text
T_b = patch_T(p)
Y_b = patch_Y(:,p)
h_b = h(T_b,Y_b,p0)
rho_b = rho(T_b,Y_b,p0)
```

Non-Dirichlet species boundaries use the interior composition for the boundary
thermodynamic state.

---

## 8. `&species_input`

| Option | Type | Description |
|---|---:|---|
| `enable_species` | logical | Master switch for species transport. |
| `enable_reactions` | logical | Reaction source switch. Current supported validation path keeps this `.false.`. |
| `enable_cantera` | logical | Use Cantera species diffusivities `D_k`. Separate from fluid Cantera switch. |
| `nspecies` | integer | Number of transported species. Must be within compiled limit. |
| `species_name(max_species)` | string array | Transported species names. Must match Cantera mechanism when Cantera is used. |
| `initial_Y(max_species)` | real array | Initial species mass fractions. |
| `species_diffusivity(max_species)` | real array | Constant fallback species diffusivities `[m^2/s]`. |

Rules:

```text
enable_reactions = .false. for current supported non-reacting runs
species mass fractions should be nonnegative
fixed-value boundary patch_Y should sum to approximately 1
```

Cantera species diffusivity:

```text
if species_input enable_cantera = .false.:
  D_k = species_diffusivity(k)

if species_input enable_cantera = .true.:
  D_k = Cantera D_k(T,Y,p0)
```

Temperature source for Cantera species transport:

```text
energy enabled:   energy%T
energy disabled:  background_temp
```

---

## 9. `&energy_input`

| Option | Type | Description |
|---|---:|---|
| `enable_energy` | logical | Master switch for sensible-enthalpy transport. |
| `enable_cantera_thermo` | logical | Use Cantera for `h <-> T`, `cp`, `lambda`, and `rho_thermo`. |
| `enable_species_enthalpy_diffusion` | logical | Include optional `-div(sum_k h_k J_k)` correction. Requires species and Cantera thermo for meaningful `h_k`. |
| `thermo_update_interval` | integer | Current supported value is `1`. |
| `thermo_default_species` | string | Bath/default species when species transport is disabled. |
| `initial_T` | real | Initial temperature `[K]`. Must be positive. |
| `energy_reference_T` | real | Sensible enthalpy reference temperature `[K]`. Must be positive. |
| `energy_reference_h` | real | Constant-cp reference enthalpy `[J/kg]`. |
| `energy_cp` | real | Constant heat capacity fallback `[J/kg/K]`. |
| `energy_lambda` | real | Constant thermal conductivity fallback `[W/m/K]`. |

### Energy convention

```text
transported state = h_sensible [J/kg]
dependent state   = T [K]
recovery          = T(h,Y,p0)
p0                = background_press
```

When species composition changes:

```text
preserve h
recover T = T(h,Y_new,p0)
```

Do not rebuild:

```text
h = h(T_old,Y_new,p0)
```

### Constant-cp mode

```text
h = energy_reference_h + energy_cp * (T - energy_reference_T)
T = energy_reference_T + (h - energy_reference_h) / energy_cp
cp = energy_cp
lambda = energy_lambda
```

### Cantera thermo mode

```text
h_sens = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
T = Cantera HP inverse from h_sens + h_abs(T_ref,Y,p0)
cp = cp(T,Y,p0)
lambda = thermal_conductivity(T,Y,p0)
rho_thermo = rho(T,Y,p0)
```

Preferred energy-step call:

```text
(T, cp, lambda, rho_thermo) = sync(h,Y,p0)
```

The sync must preserve `h`.

### Variable-density energy

In variable-density mode:

```text
energy update = conservative rho*h branch
transport%rho <- energy%rho_thermo after thermo sync
```

Primary energy diagnostics include:

```text
relative_last_energy_update_balance_defect
rel_operator_recon_defect
rel_output_recon_defect
```

---

## 10. `&output_input`

| Option | Type | Description |
|---|---:|---|
| `output_dir` | string | Root output directory. Must be non-empty. |
| `write_vtu` | logical | Enables VTU/PVTU/PVD visualization output. |
| `write_diagnostics` | logical | Enables CSV diagnostics. |

Output layout:

```text
<output_dir>/VTK/
  ParaView VTU/PVTU/PVD files

<output_dir>/diagnostics/
  CSV diagnostics
```

Common diagnostics:

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
```

Not all files appear in all modes.

---

## 11. `&profiling_input`

| Option | Type | Description |
|---|---:|---|
| `enable_profiling` | logical | Enable MPI-aware profiler. |
| `nested_profiling` | logical | Show nested profiler tree. Timings are inclusive. |
| `write_cantera_cache_stats` | logical | Print final Cantera cache call/hit/miss statistics. |
| `variable_density_debug` | logical | Enable verbose variable-density debug prints. Use only while debugging. |

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

---

## 12. Cantera demo mechanism quick table

Curated demo mechanisms under `mechanisms/` use a phase named `gas`:

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

For converted CHEMKIN mechanisms:

```bash
ck2yaml --input chem.inp --thermo therm.dat --transport tran.dat --output mechanism.yaml
```

Then select the generated phase in `case.nml`.

---

## 13. Minimal hydrodynamic pattern

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
  enable_species_enthalpy_diffusion = .false.
  thermo_update_interval = 1
/
```

Expected:

```text
constant rho
constant nu
div(u)=0
no species or energy solve
```

---

## 14. Guarded variable-density pattern

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

Expected:

```text
active density = Cantera rho_thermo
projection target = div(u)=S_projection
species transport = rho*Y branch
energy transport = rho*h branch
reactions off
```

Validate with:

```text
divu_minus_S_projection_*
integral_drho_dt_plus_div_mass_flux_dV
relative_conservative_residual_l2
relative_last_energy_update_balance_defect
rel_operator_recon_defect
rel_output_recon_defect
```

---

## 15. Common mistakes

```text
1. Treating old docs that say variable density is unsupported as current.
2. Enabling variable density without Cantera thermo energy.
3. Confusing fluid_input enable_cantera with species_input enable_cantera.
4. Confusing background_press with projection pressure.
5. Setting patch_T without patch_temperature_type="fixed_value".
6. Fixed species patch_Y not summing to approximately 1.
7. Using interior Y for fixed-temperature fixed-composition inlet enthalpy.
8. Rebuilding h from T_old after species transport instead of preserving h.
9. Expecting equal velocity to mean equal mass flux in variable-density mode.
10. Judging variable-density projection by raw max_div instead of divu_minus_S_projection_*.
11. Reading one rank .vtu file from an MPI run instead of .pvtu/.pvd.
12. Reusing old output after VTK schema changes.
```

---

## 16. Maintenance note

This file should remain a compact reference, not a patch-history log.

Put long development history in changelogs or development notes.  This reference
should answer:

```text
What option exists?
What type is it?
What values are valid?
What does it do now?
What combinations are supported?
```