title: Output Directory Layout

# Output Directory Layout

This document defines the current output-directory contract for LowMachReact-Hex.

The solver separates visualization files from diagnostic and validation data:

```text
<output_dir>/
  VTK/
  diagnostics/
```

`output_dir` is configured in `case.nml`:

```fortran
&output_input
  output_dir = "cases/example/output"
  write_vtu = .true.
  write_diagnostics = .true.
/
```

This document is intentionally current-state oriented.  It should not contain a
long patch history.  Patch notes, validation explanations, and diagnostic
interpretation belong in the relevant validation documents.

---

## 1. Top-level layout

The main run directory is:

```text
<output_dir>/
```

The solver creates the normal output subdirectories during startup:

```text
<output_dir>/VTK/
<output_dir>/diagnostics/
```

Meaning:

| Path | Purpose |
|---|---|
| `<output_dir>/VTK/` | ParaView visualization files. |
| `<output_dir>/diagnostics/` | CSV diagnostics, validation metrics, and audit-friendly scalar outputs. |

Postprocessing scripts should not assume output files live directly inside
`output_dir`.  Use:

```text
vtk_dir = output_dir / "VTK"
diagnostics_dir = output_dir / "diagnostics"
```

---

## 2. ParaView visualization files

Visualization files are written under:

```text
<output_dir>/VTK/
```

Typical files:

```text
*.vtu
*.pvtu
*.pvd
```

Meaning:

| File type | Meaning |
|---|---|
| `.vtu` | VTK XML unstructured-grid piece files. |
| `.pvtu` | Parallel VTK master file that references `.vtu` pieces. |
| `.pvd` | ParaView time-series collection file. |

Recommended ParaView workflow:

```text
Open the .pvd file for a time series.
Open the .pvtu file for a single parallel timestep.
Do not use one rank's .vtu file as the full-domain result.
```

For MPI runs, `.pvd` or `.pvtu` should be treated as the full-domain
visualization entry point.

---

## 3. Cell-centered finite-volume fields

Most solver fields are finite-volume cell-centered data.

Common visualization fields may include:

```text
velocity
pressure
rho
nu
temperature
enthalpy
qrad
cp
thermal_conductivity
rho_thermo
thermal_diffusivity
thermo_pressure
Y_<species>
D_<species>
sum_Y
mass_flux_vector
mass_flux_divergence
```

Important interpretation:

| Field | Meaning |
|---|---|
| `pressure` | Projection pressure-like field, not Cantera thermodynamic pressure. |
| `thermo_pressure` | Uniform thermodynamic/background pressure passed to Cantera. |
| `rho` | Active flow density `transport%rho`. |
| `rho_thermo` | Cantera thermodynamic density. Diagnostic in constant-density mode; active density source in variable-density mode. |
| `enthalpy` | Transported sensible enthalpy `h`. |
| `qrad` | Volumetric energy source; positive adds energy to the gas. |
| `mass_flux_vector` | Cell-data visualization of mass flux density, approximately `rho*u`. |
| `mass_flux_divergence` | Cell-data reconstruction of `div(rho u)` from face mass fluxes. |

ParaView filters such as `Plot Over Line`, `Resample To Line`, and `Cell Data to
Point Data` may interpolate or average fields.  Quantitative scripts should
compare like with like:

```text
raw CellData extraction <-> raw cell-centered finite-volume output
ParaView line sample    <-> equivalent scripted interpolation
PVTU/PVD dataset        <-> full-domain data, not one rank piece
```

---

## 4. Variable-density visualization fields

When `enable_variable_density = .true.`, additional low-Mach debug fields may be
written under `VTK/`.

Common variable-density debug fields include:

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

Interpretation:

| Field | Meaning |
|---|---|
| `lowmach_source_projection` | Source copied immediately before projection RHS assembly. |
| `lowmach_source_current` | Source after the latest energy/thermo density update. |
| `divu_minus_S_projection` | Projection-time low-Mach residual. Primary spatial projection debug field. |
| `divu_minus_S_current` | Source-evolution diagnostic, not the primary projection pass/fail metric. |
| `rho_current_minus_projection` | Difference between current density and projection-time density snapshot. |
| `continuity_residual_estimate` | Spatial diagnostic for conservative-continuity residual localization. |
| `rho_h_*` fields | Energy density/time-level reconciliation visualization fields. |

Use these fields to locate problems in space.  Use CSV diagnostics for global
validation decisions.

---

## 5. Diagnostics files

CSV diagnostics are written under:

```text
<output_dir>/diagnostics/
```

Not every file appears in every run.  Files are gated by runtime options such as:

```text
write_diagnostics
enable_energy
enable_species
enable_variable_density
enable_cantera_thermo
write_cantera_cache_stats
```

### 5.1 Core flow diagnostics

Common solver-generated files:

```text
diagnostics.csv
```

Typical contents include:

```text
step
time
dt
CFL
pressure iterations
pressure residual
max/rms divergence or low-Mach runtime residuals
net boundary flux
kinetic energy
maximum velocity
total mass
```

In constant-density mode, raw `max_div` and `rms_div` are meaningful projection
diagnostics because the target is `div(u)=0`.

In variable-density mode, raw divergence is not expected to vanish.  The target
is `div(u)=S_projection`.

### 5.2 Energy diagnostics

Energy-related solver-generated files may include:

```text
energy_diagnostics.csv
enthalpy_energy_budget.csv
```

Typical energy fields include:

```text
h_min/h_max/h_mean
T_min/T_max/T_mean
qrad_min/qrad_max
qrad_integral
rho_h_integral
reported_conductive_boundary_flux_out
relative_last_energy_update_balance_defect
rel_operator_recon_defect
rel_output_recon_defect
```

Primary energy closure metrics are:

```text
relative_last_energy_update_balance_defect
rel_operator_recon_defect
rel_output_recon_defect
```

The energy budget distinguishes direct operator closure, operator-consistent
`rho*h` budget, and output-state density reconciliation.

### 5.3 Species diagnostics

Species-related solver-generated files may include:

```text
species_energy_conservation.csv
species_integrals.csv
```

Typical species diagnostics include:

```text
transported_species_mass_sum
sum_Y_min
sum_Y_max
sum_Y_mean
sum_Y_l2
per-species integrals
per-species approximate boundary fluxes
```

The species boundary flux diagnostics are validation/trend diagnostics.  They
should be interpreted with the sign convention documented in the diagnostic
guide and with awareness of the boundary-state assumptions.

### 5.4 Variable-density diagnostics

Variable-density solver-generated files may include:

```text
variable_density_diagnostics.csv
variable_density_compatibility.csv
variable_density_transport_conservation.csv
variable_density_continuity_residual.csv
variable_density_worst_cell.csv
variable_density_worst_cell_faces.csv
variable_density_boundary_residual_summary.csv
variable_density_boundary_residual_cells_rank<RANK>.csv
variable_density_projection_audit_cells_rank<RANK>.csv
variable_density_projection_audit_faces_rank<RANK>.csv
```

Primary variable-density projection metrics:

```text
divu_minus_S_projection_max
divu_minus_S_projection_l2
relative_divu_minus_S_projection_max
relative_divu_minus_S_projection_l2
net_boundary_volume_flux_minus_integral_S_projection_dV
```

Current-source columns such as:

```text
divu_minus_S_current_*
S_current_minus_S_projection_*
net_boundary_volume_flux_minus_integral_S_current_dV
```

are source-evolution diagnostics.  They are useful, but they are not the primary
projection pass/fail metric.

Primary conservative-continuity metrics:

```text
integral_drho_dt_plus_div_mass_flux_dV
conservative_residual_l2
relative_conservative_residual_l2
mass_balance_defect_*
```

### 5.5 Cantera cache diagnostics

When Cantera cache statistics are enabled, diagnostics may include:

```text
cantera_cache_stats.csv
```

or equivalent final cache-statistics reporting, depending on the current build
and output path.

The cache diagnostics are performance/debugging data.  They are not physics
validation metrics.

---

## 6. Solver-generated versus offline-generated files

Most CSV files in `<output_dir>/diagnostics/` are written directly by the solver.

Some validation files are produced by offline tools and are not solver state
files.

### 6.1 Coupled transport audit

Offline tool:

```bash
python tools/diagnostics/audit_coupled_transport_conservation.py \
  --output-dir <output_dir> \
  --write-csv
```

Possible output:

```text
<output_dir>/diagnostics/coupled_transport_audit.csv
```

This file cross-checks existing mass, species, continuity, and reconciled energy
diagnostics.  It is derived validation output.

### 6.2 Variable-density validation checker

Offline tool:

```bash
python tools/diagnostics/check_variable_density_validation.py \
  --output-dir <output_dir>
```

or:

```bash
python tools/diagnostics/check_variable_density_validation.py \
  --diagnostics-dir <output_dir>/diagnostics
```

This checker consumes solver-generated CSV files and reports pass/fail status for
accepted validation metrics.

### 6.3 Validation matrix summary

Offline tool:

```text
tools/diagnostics/run_variable_density_validation_matrix.py
```

Possible output:

```text
validation_matrix_summary.csv
```

This file summarizes validation results across multiple output directories or
case variants.  It is not produced by a single solver run unless the offline
matrix tool is invoked.

---

## 7. Output files by mode

### 7.1 Flow-only constant-density run

Likely output:

```text
<output_dir>/VTK/*.vtu
<output_dir>/VTK/*.pvtu
<output_dir>/VTK/*.pvd
<output_dir>/diagnostics/diagnostics.csv
```

### 7.2 Energy-enabled run

Additional likely output:

```text
<output_dir>/diagnostics/energy_diagnostics.csv
<output_dir>/diagnostics/enthalpy_energy_budget.csv
```

and additional VTK fields such as:

```text
temperature
enthalpy
qrad
cp
thermal_conductivity
rho_thermo
```

### 7.3 Species-enabled run

Additional likely output:

```text
<output_dir>/diagnostics/species_energy_conservation.csv
<output_dir>/diagnostics/species_integrals.csv
```

and additional VTK fields such as:

```text
Y_<species>
D_<species>
sum_Y
```

### 7.4 Variable-density low-Mach run

Additional likely output:

```text
<output_dir>/diagnostics/variable_density_diagnostics.csv
<output_dir>/diagnostics/variable_density_compatibility.csv
<output_dir>/diagnostics/variable_density_transport_conservation.csv
<output_dir>/diagnostics/variable_density_continuity_residual.csv
```

Debug/localization files may also appear:

```text
variable_density_worst_cell.csv
variable_density_worst_cell_faces.csv
variable_density_boundary_residual_summary.csv
variable_density_boundary_residual_cells_rank<RANK>.csv
variable_density_projection_audit_cells_rank<RANK>.csv
variable_density_projection_audit_faces_rank<RANK>.csv
```

Additional VTK debug fields may appear when the variable-density output path is
enabled.

---

## 8. Sign conventions used by diagnostics

Boundary flux convention:

```text
positive boundary flux = outward from the domain
negative boundary flux = inflow to the domain
```

Mass conservation diagnostic convention:

```text
dM/dt + net_boundary_mass_flux ~= 0
```

Energy budget convention:

```text
boundary advective rho*h flux is positive outward
conductive boundary flux is positive when energy leaves the domain
qrad > 0 adds energy to the gas
qrad < 0 removes energy from the gas
```

Variable-density projection convention:

```text
constant-density mode:
  target = div(u) = 0

variable-density mode:
  target = div(u) = S_projection
```

---

## 9. Cleaning old output

After changing VTK arrays, PVTU schema, diagnostics columns, or output layout,
delete the previous output directory before rerunning:

```bash
rm -rf <output_dir>
```

Example:

```bash
rm -rf cases/counterflow_nonreacting_gmsh/output
```

Old VTK/PVTU schemas can confuse ParaView, and old CSV files can confuse
postprocessing scripts.

---

## 10. Recommended postprocessing rules

Use these rules for robust scripts:

```text
1. Treat output_dir as a run container, not as a flat file directory.
2. Read visualization files from output_dir/VTK.
3. Read scalar diagnostics from output_dir/diagnostics.
4. Use .pvd or .pvtu for full-domain visualization.
5. Do not use one rank .vtu file as the full MPI result.
6. Do not assume every diagnostics file exists in every mode.
7. Check runtime options before expecting variable-density, energy, species, or
   Cantera diagnostics.
8. Prefer CSV diagnostics for global conservation and closure decisions.
9. Use VTK fields for spatial localization and debugging.
10. Delete old output after output-schema changes.
```

Python path pattern:

```python
from pathlib import Path

output_dir = Path("cases/example/output")
vtk_dir = output_dir / "VTK"
diagnostics_dir = output_dir / "diagnostics"
```

---

## 11. Related documentation

See also:

```text
doc_src/paraview_output_fields.md
doc_src/validation_metrics.md
doc_src/validation_automation.md
doc_src/input_configuration_guide.md
doc_src/case_nml_reference.md
doc_src/numerical_method.md
```

Use:

```text
paraview_output_fields.md
```

for the meaning of individual VTK arrays.

Use:

```text
validation_metrics.md
```

to distinguish primary pass/fail metrics from explanatory or reconciliation
columns.

Use:

```text
validation_automation.md
```

for offline checker and validation-matrix workflows.

---

## 12. Maintenance rule

Keep this document focused on the output layout and file meaning.

Do not append long patch histories here.  This document should answer:

```text
Where are files written?
Which files are solver-generated?
Which files are offline-derived?
Which files appear only in certain modes?
Which files should users open in ParaView?
Which CSV files should scripts read?
```