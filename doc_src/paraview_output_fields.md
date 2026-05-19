title: ParaView Output Field Guide


# ParaView Output Field Guide

This guide explains the main quantities written to VTU/PVTU output for viewing
and debugging LowMachReact-Hex results in ParaView.

ParaView fields are primarily **spatial localization tools**.  For global
conservation, projection quality, and pass/fail validation, use the CSV
diagnostics under:

```text
<output_dir>/diagnostics/
```

The most important distinction is:

```text
constant-density mode:
  projection target is div(u) = 0

variable-density low-Mach mode:
  projection target is div(u) = S_projection
```

Therefore, in variable-density mode, raw `div(u)` is not by itself a projection
error.  The projection-error field is:

```text
divu_minus_S_projection = div(u) - S_projection
```

---

## 1. Output location

ParaView files are written under:

```text
<output_dir>/VTK/
```

Typical files:

```text
*.vtu
*.pvtu
*.pvd
```

Recommended usage:

```text
Open *.pvd for a time series.
Open *.pvtu for a full parallel timestep.
Do not interpret one rank's *.vtu file as the full MPI result.
```

Most arrays are finite-volume **cell-centered** fields.  ParaView operations such
as `Plot Over Line`, `Resample To Line`, and `Cell Data to Point Data` may
interpolate or average them.

For quantitative comparison, compare like with like:

```text
raw CellData extraction <-> raw finite-volume cell data
ParaView line sample    <-> equivalent scripted interpolation
PVTU/PVD dataset        <-> full-domain result, not one rank piece
```

---

## 2. CSV-first validation rule

Use ParaView to answer:

```text
Where is the residual large?
Where is the density changing?
Where is the source localized?
Where is the energy reconciliation localized?
```

Use CSV diagnostics to answer:

```text
Did the projection converge?
Does conservative continuity close?
Does the energy update close?
Does the reconciled rho*h budget close?
```

Primary variable-density validation metrics are documented in
`doc_src/validation_metrics.md`.

Common examples:

| Category | Primary CSV metrics |
|---|---|
| Projection | `divu_minus_S_projection_*` |
| Conservative continuity | `integral_drho_dt_plus_div_mass_flux_dV`, `relative_conservative_residual_l2` |
| Direct energy update | `relative_last_energy_update_balance_defect` |
| Reconciled rho*h budget | `rel_operator_recon_defect`, `rel_output_recon_defect` |

Current-source residuals, raw `div(u)`, and unreconciled `rho*h` budget columns
are explanatory diagnostics.  They should not replace the primary metrics above.

---

## 3. Core flow fields

Common flow fields may include:

| Field | Meaning | Units | Notes |
|---|---|---:|---|
| `velocity` or velocity components | Cell-centered velocity. | m/s | Primary flow variable. |
| `pressure` | Projection pressure-like field. | solver pressure units | Not the Cantera thermodynamic pressure. |
| `rho` | Active flow density `transport%rho`. | kg/m^3 | Constant in constant-density mode; Cantera-synced in variable-density mode. |
| `rho_current` | Current active density. | kg/m^3 | Usually available in variable-density debug output. |
| `nu` | Kinematic viscosity. | m^2/s | Constant unless variable-viscosity mode is active. |
| `mass_flux_vector` | Cell-centered visualization of approximately `rho*u`. | kg/(m^2 s) | The exact conservative mass flux used by the solver is face-centered. |
| `mass_flux_divergence` | Reconstructed `div(rho u)` from face mass fluxes. | kg/(m^3 s) | Useful for localizing mass-transport behavior. |

Important pressure distinction:

```text
pressure:
  hydrodynamic/projection pressure-like field

background_press / p0:
  thermodynamic pressure passed to Cantera
```

Do not interpret `pressure` as the Cantera EOS pressure.

---

## 4. Constant-density projection fields

For constant-density incompressible/projection runs:

```text
target: div(u) = 0
```

Useful checks:

```text
max |div(u)|
rms div(u)
net boundary volume flux
```

If a divergence visualization field is written, large localized values usually
indicate a projection, boundary, or flux-balance issue.

---

## 5. Variable-density low-Mach source fields

For variable-density low-Mach runs:

```text
target: div(u) = S_projection
```

Raw divergence measures volumetric expansion or contraction and is allowed to be
nonzero.

Common fields:

| Field | Meaning | How to interpret |
|---|---|---|
| `lowmach_source_projection` | Source copied immediately before pressure RHS assembly. | Correct source target for the projected velocity field. |
| `lowmach_source_current` | Source after the latest energy/thermo density sync. | Useful for visualizing source evolution after projection. |
| `lowmach_source_difference` | `S_current - S_projection`. | Time-level/source-evolution difference. |
| `divu_recomputed` | Recomputed `div(u)` from corrected volumetric face fluxes. | Physical expansion/contraction; not an error by itself in variable-density mode. |
| `divu_minus_S_projection` | `div(u) - S_projection`. | Primary ParaView projection-error field in variable-density mode. |
| `divu_minus_S_current` | `div(u) - S_current`. | Source-evolution diagnostic after energy/thermo has advanced the source. |

Recommended projection-debug field:

```text
divu_minus_S_projection
```

Recommended source-time-level fields:

```text
lowmach_source_difference
divu_minus_S_current
```

Interpretation:

```text
large divu_recomputed + small divu_minus_S_projection:
  correct variable-density expansion, not a projection failure

large divu_minus_S_projection:
  possible projection, pressure solve, boundary flux, or compatibility issue

large divu_minus_S_current + small divu_minus_S_projection:
  source changed after projection; expected time-level drift
```

---

## 6. Conservative low-Mach source decomposition

The active variable-density low-Mach source is:

```text
S = (rho_old - rho) / (rho * dt) - (u.grad(rho)) / rho
```

This is often interpreted as:

```text
S = S_history + S_advective_density
```

where:

```text
S_history           = (rho_old - rho) / (rho * dt)
S_advective_density = -(u.grad(rho)) / rho
```

Common fields:

| Field | Meaning | How to interpret |
|---|---|---|
| `lowmach_source_history_estimate` | Estimate of `(rho_old-rho)/(rho*dt)`. | Local expansion/contraction from density history alone. |
| `u_dot_grad_rho` | Explicit density-advection term `u.grad(rho)`. | Shows where velocity carries density gradients. |
| `lowmach_source_advective_density` | `-(u.grad(rho))/rho`. | Advective-density contribution to the low-Mach source. |
| `continuity_residual_estimate` | Estimated local conservative-continuity residual. | Spatial localization only; use CSV for global closure. |

Useful identity:

```text
div(rho u) = rho div(u) + u.grad(rho)
```

This decomposition is useful for distinguishing two causes of local expansion:

```text
density history:
  local thermodynamic density changed between time levels

density advection:
  flow advected density gradients through the cell
```

---

## 7. Density time-level fields

Variable-density output may include density snapshots from different time levels.

| Field | Meaning | How to interpret |
|---|---|---|
| `rho_current` | Current active density after thermo/EOS sync. | Density used by current transport/conservation diagnostics. |
| `rho_projection` | Density snapshot associated with the projection/source time level. | Compare with `rho_current` to see post-projection density evolution. |
| `rho_current_minus_projection` | `rho_current - rho_projection`. | Density time-level drift. |

Useful fields to plot together:

```text
rho_current
rho_projection
rho_current_minus_projection
lowmach_source_difference
divu_minus_S_current
```

If `rho_current_minus_projection` is large in a region, `lowmach_source_difference`
and `divu_minus_S_current` are often large there as well.

---

## 8. Conservative mass-flux fields

The exact finite-volume mass flux is face-centered.  The VTK writer exposes
cell-centered reconstructions for visualization.

| Field | Meaning | How to interpret |
|---|---|---|
| `mass_flux_vector` | Cell-centered visualization of `rho*u`. | Useful for mass-flow direction and magnitude; not the exact face flux array. |
| `mass_flux_divergence` | Reconstructed `div(rho u)`. | General density-weighted flux-divergence visualization. |
| `mass_flux_divergence_recomputed` | Recomputed `div(rho u)` from face mass fluxes. | Useful with continuity residual fields. |

Use CSV diagnostics for global mass validation:

```text
variable_density_transport_conservation.csv
variable_density_continuity_residual.csv
```

---

## 9. Energy and thermodynamic fields

Common energy/thermo fields may include:

| Field | Meaning | Units | Notes |
|---|---|---:|---|
| `temperature` or `T` | Gas temperature. | K | Recovered from `h,Y,p0` in Cantera thermo mode. |
| `enthalpy` or `h` | Transported sensible enthalpy. | J/kg | Primary energy state. |
| `cp` | Mixture heat capacity. | J/(kg K) | From Cantera or constant model. |
| `thermal_conductivity` or `lambda` | Thermal conductivity. | W/(m K) | From Cantera or constant model. |
| `thermal_diffusivity` | Thermal diffusivity-style derived field, when written. | m^2/s | Diagnostic/visualization field. |
| `rho_thermo` | Thermodynamic density from Cantera/EOS. | kg/m^3 | Diagnostic in constant-density mode; active density source in variable-density mode. |
| `qrad` | Volumetric source storage. | W/m^3 | Positive adds energy to gas. |
| `species_enthalpy_diffusion` | Species-enthalpy diffusion contribution, when written. | implementation-dependent diagnostic units | Present only when the corresponding diagnostic/output path is enabled. |

### 9.1 Thermodynamic pressure field

If a thermodynamic-pressure field such as `thermo_pressure` is written, interpret
it as:

```text
thermo_pressure = background_press = p0
```

It is the uniform Cantera/EOS pressure used by calls such as:

```text
rho(T,Y,p0)
T(h,Y,p0)
cp(T,Y,p0)
lambda(T,Y,p0)
```

It is not the projection pressure field.

In the current low-Mach formulation, `p0` is spatially uniform.  A repeated
per-cell `thermo_pressure` field is therefore metadata-like; it is useful mainly
to avoid confusing projection pressure with thermodynamic pressure.

---

## 10. Energy-density reconciliation fields

Variable-density energy output may include local fields that explain the
difference between the output-state `rho*h` snapshot and the density/time level
used internally by the energy operator.

| Field | Meaning | Units | Notes |
|---|---|---:|---|
| `rho_h_output_state` | Output-snapshot energy density, approximately `transport%rho * energy%h`. | J/m^3 | Local state represented by the output fields. |
| `rho_h_operator_consistent` | Cellwise `rho*h` energy density at the energy-operator density/time level. | J/m^3 | Populated during the energy update; before the first update may fall back to output state. |
| `rho_h_density_reconciliation` | `rho_h_output_state - rho_h_operator_consistent`. | J/m^3 | Signed local bookkeeping/time-level difference. |
| `relative_rho_h_density_reconciliation` | Relative local reconciliation magnitude. | nondim | Helps locate where the reconciliation is large relative to local energy density. |

Recommended fields:

```text
rho_h_density_reconciliation
relative_rho_h_density_reconciliation
```

Do not interpret these as standalone conservation errors.  Confirm global energy
closure in:

```text
<output_dir>/diagnostics/enthalpy_energy_budget.csv
```

Primary global energy metrics include:

```text
relative_last_energy_update_balance_defect
rel_operator_recon_defect
rel_output_recon_defect
```

---

## 11. Species fields

Species fields are written as cell-centered mass fractions, usually one field
per transported species.

Common naming patterns:

```text
Y_CH4
Y_O2
Y_N2
Y_<species>
sum_Y
D_<species>
```

Interpretation:

| Field | Meaning |
|---|---|
| `Y_<species>` | Mass fraction of transported species. |
| `sum_Y` | Local sum of transported mass fractions. |
| `D_<species>` | Species diffusivity used for that transported species, when written. |

Species fields influence Cantera thermodynamics when Cantera thermo is enabled:

```text
rho(T,Y,p0)
h(T,Y,p0)
T(h,Y,p0)
cp(T,Y,p0)
lambda(T,Y,p0)
```

In variable-density mode, species gradients can create density gradients and
therefore low-Mach expansion source terms.

---

## 12. Recommended ParaView workflows

### 12.1 Projection validation localization

Plot:

```text
divu_minus_S_projection
```

Also inspect:

```text
lowmach_source_projection
divu_recomputed
```

Expected:

```text
divu_minus_S_projection should be near solver/local flux tolerance.
```

### 12.2 Source time-level drift

Plot:

```text
lowmach_source_difference
divu_minus_S_current
rho_current_minus_projection
```

Expected:

```text
These fields should correlate when energy/thermo updates the density/source
after projection.
```

### 12.3 Conservative continuity localization

Plot:

```text
continuity_residual_estimate
mass_flux_divergence_recomputed
u_dot_grad_rho
lowmach_source_advective_density
```

Then confirm global norms and integrals with CSV diagnostics.

### 12.4 Density-expansion physics

Plot:

```text
rho_current
temperature
lowmach_source_history_estimate
lowmach_source_advective_density
lowmach_source_current
```

This helps distinguish expansion caused by local density change from expansion
caused by advection through density gradients.

### 12.5 Energy budget reconciliation

Plot:

```text
rho_h_density_reconciliation
relative_rho_h_density_reconciliation
```

Then confirm global closure in:

```text
diagnostics/enthalpy_energy_budget.csv
```

### 12.6 Species mixing

Plot:

```text
Y_<species>
sum_Y
rho_current
temperature
```

This helps see whether composition gradients are driving density, temperature,
or low-Mach source changes.

---

## 13. Quick interpretation table

| What you see | Likely meaning |
|---|---|
| large `divu_recomputed`, small `divu_minus_S_projection` | Correct variable-density expansion; not a projection error. |
| large `divu_minus_S_projection` | Projection, pressure solve, boundary flux, or compatibility issue. |
| large `divu_minus_S_current`, small `divu_minus_S_projection` | Source changed after projection; expected time-level drift. |
| large `lowmach_source_advective_density` | Density advection is important locally. |
| large `continuity_residual_estimate` | Possible local conservative mass-balance issue; confirm with CSV diagnostics. |
| large `rho_h_density_reconciliation` | Local output-state/operator-density bookkeeping difference; confirm global energy closure with CSV. |
| uniform `thermo_pressure` | Uniform `background_press`; not the projection pressure. |
| unexpected profile difference between ParaView line plot and script extraction | Likely interpolation/cell-data mismatch or single-rank VTU vs full PVTU/PVD mismatch. |

---

## 14. Common mistakes

```text
1. Treating raw div(u) as an error in variable-density mode.
2. Validating from ParaView fields instead of CSV diagnostics.
3. Opening a single-rank .vtu file from an MPI run.
4. Comparing ParaView line-sampled data against raw cell centers.
5. Confusing projection pressure with thermodynamic pressure.
6. Treating rho_thermo as active density in constant-density mode.
7. Treating rho_h_density_reconciliation as a standalone conservation error.
8. Expecting every field listed here to appear in every run.
```

---

## 15. Related files

See also:

```text
doc_src/output_layout.md
doc_src/validation_metrics.md
doc_src/validation_automation.md
doc_src/input_configuration_guide.md
doc_src/case_nml_reference.md
doc_src/numerical_method.md
```

Use `output_layout.md` for directory structure, `validation_metrics.md` for
primary pass/fail metrics, and `numerical_method.md` for the equations behind
the fields.

---

## 16. Maintenance rule

Keep this document focused on the meaning of VTK/PVTU fields.

Do not append long patch histories here.  This guide should answer:

```text
What does this ParaView field mean?
What are its units?
Is it a primary error metric or a localization aid?
Which CSV file should confirm the global result?
Which fields should be plotted together?
```