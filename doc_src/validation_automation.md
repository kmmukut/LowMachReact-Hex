title: Validation Automation Runbook


# Validation Automation Runbook

This runbook documents the offline validation automation used for the
variable-density low-Mach validation track in LowMachReact-Hex.

The validation tools do **not** change solver numerics.  They read solver
diagnostic CSV files, run accepted checks, and summarize validation status as
`PASS`, `WARN`, `FAIL`, or `SKIP`.

The tools are intended to answer:

```text
Did this output directory pass the accepted variable-density checks?
Did projection, continuity, and energy closure remain acceptable?
Did coupled mass/species/enthalpy consistency pass?
Which cases in a validation matrix need attention?
```

They are not a replacement for inspecting the underlying diagnostics when a case
warns or fails.

---

## 1. Scope

This runbook covers offline validation for the current guarded non-reacting
variable-density low-Mach path:

```text
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
thermo_update_interval = 1
enable_reactions = .false.
```

The automation is also useful for checking diagnostic completeness and energy
closure outputs, but the primary target is the variable-density validation
workflow.

The tools do not define new physics criteria.  Accepted primary metrics are
documented in:

```text
doc_src/validation_metrics.md
```

---

## 2. Tool chain

The current validation automation stack is:

```text
tools/diagnostics/check_variable_density_validation.py
tools/diagnostics/audit_coupled_transport_conservation.py
tools/diagnostics/run_variable_density_validation_matrix.py
```

Roles:

| Tool | Role |
|---|---|
| `check_variable_density_validation.py` | One-output-directory checker for accepted projection, continuity, and energy validation metrics. |
| `audit_coupled_transport_conservation.py` | Offline coupled mass/species/enthalpy consistency audit. |
| `run_variable_density_validation_matrix.py` | Matrix runner that invokes the checker and optional audit over one or more output directories or case rows. |

These are offline tools.  They consume files already written by the solver unless
the matrix runner is explicitly told to run commands.

---

## 3. Solver outputs versus offline outputs

The solver writes runtime outputs under each case output directory:

```text
<output_dir>/VTK/
<output_dir>/diagnostics/
```

Examples of solver-produced diagnostics include:

```text
<output_dir>/diagnostics/diagnostics.csv
<output_dir>/diagnostics/energy_diagnostics.csv
<output_dir>/diagnostics/enthalpy_energy_budget.csv
<output_dir>/diagnostics/variable_density_compatibility.csv
<output_dir>/diagnostics/variable_density_continuity_residual.csv
<output_dir>/diagnostics/variable_density_transport_conservation.csv
<output_dir>/diagnostics/species_energy_conservation.csv
<output_dir>/diagnostics/species_integrals.csv
```

The offline validation tools may write derived reports such as:

```text
<output_dir>/diagnostics/coupled_transport_audit.csv
validation_matrix_summary.csv
```

These derived files are validation reports.  They are not solver state and are
not used by the solver to continue a run.

---

## 4. One-case existing-output validation

Use this when a case has already been run and its diagnostic CSV files exist.

Run the accepted variable-density checker:

```bash
python tools/diagnostics/check_variable_density_validation.py \
  --output-dir cases/rectangle_2D/output
```

Run the coupled transport audit and write its CSV report:

```bash
python tools/diagnostics/audit_coupled_transport_conservation.py \
  --output-dir cases/rectangle_2D/output \
  --write-csv
```

Run both through the matrix runner using a direct `--case` entry:

```bash
python tools/diagnostics/run_variable_density_validation_matrix.py \
  --case rectangle_2D=cases/rectangle_2D/output \
  --write-audit-csv \
  --summary-csv validation_matrix_summary.csv
```

This is the recommended quick validation command for an existing single output
directory.

---

## 5. Matrix validation over existing outputs

The matrix runner can read a CSV or JSON matrix.  The example CSV is:

```text
tools/diagnostics/variable_density_validation_matrix.example.csv
```

Run the example matrix over existing outputs:

```bash
python tools/diagnostics/run_variable_density_validation_matrix.py \
  --matrix tools/diagnostics/variable_density_validation_matrix.example.csv \
  --write-audit-csv \
  --summary-csv validation_matrix_summary.csv
```

By default, `run_command` entries are not executed.  This makes the command safe
for auditing existing outputs without accidentally rerunning simulations.

---

## 6. Matrix CSV format

The matrix CSV uses these columns:

```text
name
output_dir
run_command
cwd
enabled
energy_only
```

Meaning:

| Column | Meaning |
|---|---|
| `name` | User-facing case or variant name in the summary. |
| `output_dir` | Path to the solver output directory for that case. |
| `run_command` | Optional command to run before validation when `--run-commands` is used. |
| `cwd` | Optional working directory for `run_command`. |
| `enabled` | Whether this row is included. Disabled rows are skipped. |
| `energy_only` | Marks rows intended for energy-focused checking when supported by the checker. |

Recommended workflow for adding a new case:

```text
1. Add the row with enabled=false.
2. Run the case manually.
3. Validate the case with the one-case command.
4. Inspect WARN/FAIL details if any.
5. Enable the row only after the single-case checks pass.
```

---

## 7. Optional run-command mode

To execute each enabled row's `run_command` before validation, pass
`--run-commands`:

```bash
python tools/diagnostics/run_variable_density_validation_matrix.py \
  --matrix tools/diagnostics/variable_density_validation_matrix.example.csv \
  --run-commands \
  --write-audit-csv \
  --summary-csv validation_matrix_summary.csv
```

Use this mode for explicit validation campaigns, not for quick inspection of
existing outputs.

Before using `--run-commands`, check that:

```text
- each run_command is correct
- each cwd is correct
- old output is removed if the case should start clean
- the intended MPI rank count is encoded in the command
- the active Conda/MPI environment is correct
```

---

## 8. Strict CI-style mode

Normal matrix mode returns nonzero only if an enabled case fails.  A warning does
not fail the command unless `--fail-on-warn` is used.

Strict mode:

```bash
python tools/diagnostics/run_variable_density_validation_matrix.py \
  --matrix tools/diagnostics/variable_density_validation_matrix.example.csv \
  --write-audit-csv \
  --summary-csv validation_matrix_summary.csv \
  --fail-on-warn
```

Recommended usage:

```text
developer/local validation:
  omit --fail-on-warn unless investigating tolerance drift

pre-merge validation:
  use --fail-on-warn once the matrix cases are stable

exploratory new case:
  do not use --fail-on-warn initially; inspect WARN before promoting to CI
```

---

## 9. Status and exit-code policy

The validation tools report:

| Status | Meaning |
|---|---|
| `PASS` | All enabled required checks are within accepted thresholds. |
| `WARN` | No hard failure, but at least one metric is outside the preferred pass band. |
| `FAIL` | At least one required metric is outside the accepted failure threshold, or a required diagnostic file/tool is missing. |
| `SKIP` | A matrix row was disabled or an optional check was intentionally skipped. |

Exit-code policy:

```text
0
  No enabled case failed.  WARN still returns 0 unless --fail-on-warn is used.

nonzero
  At least one enabled case failed, or --fail-on-warn was used and at least one
  enabled case warned.
```

This policy lets exploratory validation surface warnings without breaking local
workflows, while still supporting strict pre-merge behavior.

---

## 10. Accepted metric source of truth

The automation follows the accepted metric definitions in:

```text
doc_src/validation_metrics.md
```

The matrix runner orchestrates checks.  It does not define new physics criteria.

Important interpretation:

```text
Large unreconciled rho*h defects are not validation failures when the reconciled
metrics pass.
```

Primary reconciled energy metrics include:

```text
rel_output_recon_defect
rel_operator_recon_defect
relative_last_energy_update_balance_defect
```

Primary projection metrics use the projection-time source:

```text
divu_minus_S_projection_*
```

not raw divergence and not current-source residuals.

Primary conservative continuity metrics include:

```text
integral_drho_dt_plus_div_mass_flux_dV
relative_conservative_residual_l2
```

---

## 11. Reading the summary CSV

When `--summary-csv` is provided, the matrix runner writes a summary such as:

```text
validation_matrix_summary.csv
```

Use the summary to identify which case needs attention, then inspect the
case-specific diagnostics under:

```text
<output_dir>/diagnostics/
```

Do not treat the summary as the only record.  The detailed CSV files remain the
source for metric values, trends, and time histories.

Recommended inspection order after a `WARN` or `FAIL`:

```text
1. Read validation_matrix_summary.csv.
2. Open the relevant case's diagnostics directory.
3. Inspect checker output from the terminal log.
4. Inspect variable_density_compatibility.csv for projection metrics.
5. Inspect variable_density_continuity_residual.csv for conservative continuity.
6. Inspect enthalpy_energy_budget.csv for energy closure/reconciliation.
7. Inspect coupled_transport_audit.csv if the audit was written.
8. Use ParaView only after CSVs identify the failing category.
```

---

## 12. Coupled transport audit

The coupled audit is produced by:

```bash
python tools/diagnostics/audit_coupled_transport_conservation.py \
  --output-dir <output_dir> \
  --write-csv
```

The optional output file is:

```text
<output_dir>/diagnostics/coupled_transport_audit.csv
```

The audit cross-checks existing diagnostics for consistency across:

```text
mass
species
conservative continuity
enthalpy / rho*h energy budget
```

It is intended for validation and CI-style postprocessing.  It is not ParaView
output and is not used by the solver.

---

## 13. Baseline expectation

The current documented baseline case is:

```text
rectangle_2D
```

Typical one-case validation command:

```bash
python tools/diagnostics/run_variable_density_validation_matrix.py \
  --case rectangle_2D=cases/rectangle_2D/output \
  --write-audit-csv \
  --summary-csv validation_matrix_summary.csv
```

Expected stable result:

```text
rectangle_2D  PASS
```

If this baseline does not pass, do not promote new matrix variants until the
baseline issue is understood.

---

## 14. Future matrix expansion policy

Expand the matrix one case or variant at a time.

Suggested future rows:

```text
smaller dt
larger dt
different output_interval
variable_nu on/off
species-enthalpy diffusion on/off
MPI rank-count variation
hot inlet / stronger temperature contrast
stronger species contrast
different EOS phase
different background_press
different pressure/outlet configuration
mesh refinement
```

Recommended promotion path:

```text
1. Add row as disabled.
2. Run the case manually or with --run-commands in a controlled campaign.
3. Run one-case validation.
4. Inspect all WARN/FAIL output.
5. Update code, diagnostics, tolerance policy, or docs only if needed.
6. Enable the row after it is stable.
7. Include it in strict --fail-on-warn mode only after repeated clean passes.
```

No solver patch is required just to add disabled matrix rows.  A solver patch is
needed only if an enabled variant exposes a real issue that requires a numerical,
diagnostic, tolerance-policy, or documentation change.

---

## 15. Minimal local validation checklist

After running a variable-density validation case:

```bash
python tools/diagnostics/run_variable_density_validation_matrix.py \
  --case rectangle_2D=cases/rectangle_2D/output \
  --write-audit-csv \
  --summary-csv validation_matrix_summary.csv
```

Before promoting a stable validation matrix:

```bash
python tools/diagnostics/run_variable_density_validation_matrix.py \
  --matrix tools/diagnostics/variable_density_validation_matrix.example.csv \
  --write-audit-csv \
  --summary-csv validation_matrix_summary.csv \
  --fail-on-warn
```

Checklist:

```text
- summary CSV exists
- no enabled case is FAIL
- no enabled case is WARN in strict mode
- coupled_transport_audit.csv is written when requested
- detailed diagnostics are present under each output_dir/diagnostics
- terminal output is saved for review if this is a pre-merge run
```

---

## 16. Common mistakes

```text
1. Running the matrix checker before the solver has produced diagnostics.
2. Passing a case root directory instead of the output directory.
3. Forgetting that run_command is skipped unless --run-commands is used.
4. Using --fail-on-warn for exploratory cases before tolerances are understood.
5. Treating unreconciled rho*h defects as failures when reconciled metrics pass.
6. Treating raw div(u) as the projection metric in variable-density mode.
7. Forgetting to remove old output before rerunning a case with changed diagnostics.
8. Enabling a new matrix row before the one-case validation passes.
```

---

## 17. Related documentation

See also:

```text
doc_src/validation_metrics.md
doc_src/output_layout.md
doc_src/paraview_output_fields.md
doc_src/numerical_method.md
doc_src/input_configuration_guide.md
doc_src/patch_workflow_and_next_steps.md
```

Use `validation_metrics.md` for accepted thresholds and metric meanings.

Use `output_layout.md` for which files are solver-generated versus
offline-generated.

Use `paraview_output_fields.md` for spatial localization fields after a CSV
metric identifies a problem category.

---

## 18. Maintenance rule

Keep this runbook focused on validation automation.

Do not append long patch histories here.  This document should answer:

```text
Which tools are used?
What files do they read and write?
How do I validate one existing output directory?
How do I validate a matrix?
When are run commands executed?
What do PASS/WARN/FAIL/SKIP mean?
What should I inspect after a warning or failure?
```