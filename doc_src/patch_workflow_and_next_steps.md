title: Patch Workflow and Next Steps


# Patch Workflow and Next Steps

This document defines the current development workflow for small code changes,
documentation updates, diagnostics, and validation work in LowMachReact-Hex.

It replaces older patch notes that referred to early profiling and Cantera
thermo-sync patches as the immediate next step.  Those historical patch notes
should stay in changelogs or development history, not in this current workflow
document.

The current code state is:

```text
stable baseline:
  constant-density incompressible / low-Mach projection solver

active guarded path:
  non-reacting variable-density low-Mach solver using Cantera density,
  conservative rho*Y and rho*h transport, and diagnostic-rich validation

future work:
  restart, density-aware boundary states, mass-flow boundaries, broader
  validation, reactions, heat release, and radiation coupling through qrad
```

---

## 1. Patch philosophy

Use one small patch per logical change.

A good patch should answer:

```text
What changed?
Why was it needed?
Which code paths are affected?
Which code paths should remain unchanged?
How can the change be validated?
How can it be rolled back?
```

Prefer a sequence of narrow patches over one large patch that changes numerics,
diagnostics, output, and documentation at the same time.

---

## 2. Python patch script requirements

Continue using Python patch scripts for mechanical edits.

Each patch script should:

```text
- support --dry-run and --apply
- default to dry-run behavior unless --apply is explicitly passed
- create timestamped backups under .backups/YYYYMMDD_HHMMSS/
- fail loudly when anchors are missing
- fail loudly when anchors are ambiguous
- be idempotent where practical
- print a clear summary of changed files
- print a clear summary of skipped files
- list required build and validation next steps
- avoid silently rewriting unrelated regions
```

Recommended command pattern:

```bash
python patches/patch_XXX_short_name.py --dry-run
python patches/patch_XXX_short_name.py --apply
```

For patches that modify many files, include a final summary like:

```text
changed:
  src/mod_energy.f90
  src/mod_output.f90

skipped:
  doc_src/numerical_method.md already contains target text

next:
  make BUILD=debug
  make BUILD=release
  run <validation cases>
```

---

## 3. Patch safety rules

### 3.1 Anchor discipline

Patch scripts should use precise anchors.

Good anchors:

```text
- nearby procedure name
- nearby comment plus unique code line
- complete old block with enough context
- exact public interface block
```

Bad anchors:

```text
- a single common variable name
- a repeated comment
- loose regex that can match multiple routines
- line-number-only edits
```

If an anchor matches zero or more than one location, stop and print the problem.

### 3.2 Idempotence

Where practical, a patch should detect whether the desired change is already
present and skip cleanly.

Example behavior:

```text
src/mod_output.f90: already patched, skipped
src/mod_energy.f90: applied
```

Do not duplicate blocks on repeated `--apply`.

### 3.3 Backups

Before applying edits, create backups under:

```text
.backups/YYYYMMDD_HHMMSS/
```

Keep the backup path printed in the patch summary.

Backup only files that will be changed.  Do not back up the entire repository
unless the patch truly requires it.

### 3.4 Dry-run output

`--dry-run` should report what would change without modifying files.

Dry-run should still perform:

```text
- file existence checks
- anchor checks
- ambiguity checks
- planned change summary
```

---

## 4. Classification of patches

Every patch should be classified before implementation.

### 4.1 Documentation-only patch

Examples:

```text
- update numerical_method.md
- update case_nml_reference.md
- fix FORD comments
- clarify ParaView field meanings
```

Validation:

```text
- docs build or render check
- no code build required unless documentation generation is part of CI
```

### 4.2 Diagnostics-only patch

Examples:

```text
- add CSV columns
- add ParaView fields
- add profiling timers
- add offline checker output
```

Validation:

```text
- build debug and release
- run at least one case that exercises the diagnostics
- check CSV headers and first/final rows
- check PVTU schema if VTK arrays are added
```

Diagnostics-only patches should not change numerical results except for tiny
roundoff-level effects from additional reductions or instrumentation.

### 4.3 Numerical patch

Examples:

```text
- change projection RHS
- change mass flux formation
- change species update
- change energy update
- change low-Mach source
- change boundary flux treatment
```

Validation:

```text
- build debug and release
- run baseline constant-density regression
- run affected physics cases
- compare primary diagnostics before/after
- document the expected metric movement
```

Numerical patches need the most conservative review.

### 4.4 Infrastructure patch

Examples:

```text
- build system update
- environment update
- restart I/O
- mesh converter change
- launcher update
```

Validation:

```text
- clean build
- smoke test
- MPI run
- relevant tool test
```

### 4.5 Validation-tool patch

Examples:

```text
- update check_variable_density_validation.py
- update audit_coupled_transport_conservation.py
- update validation matrix runner
```

Validation:

```text
- run tool on known passing output
- run tool on intentionally incomplete/missing diagnostics if supported
- check exit codes
- check generated CSV/report columns
```

---

## 5. Required pre-patch checklist

Before writing a patch, identify:

```text
1. Target files.
2. Affected runtime modes.
3. Whether the change is documentation-only, diagnostics-only, numerical,
   infrastructure, or validation-tooling.
4. Whether constant-density behavior should remain bitwise/diagnostically
   unchanged.
5. Whether variable-density behavior is affected.
6. Which diagnostics should move, and in what direction.
7. Which validation cases must be run.
8. Which documentation must be updated.
```

Affected mode labels:

```text
constant-density flow
constant-density species
constant-density energy
Cantera-assisted constant-density thermo
variable-density low-Mach
variable-density species rho*Y
variable-density energy rho*h
Cantera transport
Cantera thermo
VTK/PVTU output
CSV diagnostics
offline validation tools
```

---

## 6. Required post-patch checklist

After applying a patch:

```bash
git diff
```

Then check:

```text
- Did only intended files change?
- Did the patch duplicate any blocks?
- Did it preserve constant-density behavior when expected?
- Did it preserve Option A enthalpy/species coupling?
- Did it preserve background_press versus projection pressure semantics?
- Did it preserve projection-source versus current-source diagnostics?
- Did it update documentation if user-facing behavior changed?
```

Build:

```bash
make BUILD=debug
make BUILD=release
```

Run at least the validation subset appropriate to the patch type.

For output-schema changes, remove old output before rerunning:

```bash
rm -rf <output_dir>
```

---

## 7. Minimum validation ladder

Use the smallest validation ladder that actually exercises the patch.

### 7.1 Baseline flow

```text
enable_energy = .false.
enable_species = .false.
enable_variable_density = .false.
```

Purpose:

```text
protect the stable constant-density projection baseline
```

### 7.2 Constant-density species

```text
enable_species = .true.
enable_reactions = .false.
enable_variable_density = .false.
```

Check:

```text
sum_Y
boundedness
species integrals
no-flux/periodic conservation where applicable
```

### 7.3 Constant-property energy

```text
enable_energy = .true.
enable_cantera_thermo = .false.
```

Check:

```text
h/T constant-cp relation
conduction behavior
qrad sign if source is used
```

### 7.4 Cantera thermo without species

```text
enable_energy = .true.
enable_cantera_thermo = .true.
enable_species = .false.
```

Check:

```text
thermo_default_species behavior
T-h-T recovery
cp/lambda/rho_thermo output
```

### 7.5 Cantera thermo with species

```text
enable_species = .true.
enable_energy = .true.
enable_cantera_thermo = .true.
enable_reactions = .false.
```

Check:

```text
Option A: preserve h, recover T(h,Y,p0)
boundary h(T_b,Y_b,p0) behavior
species-enthalpy diffusion on/off if affected
```

### 7.6 Guarded variable-density low-Mach

```text
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
thermo_update_interval = 1
enable_reactions = .false.
```

Primary checks:

```text
divu_minus_S_projection_*
integral_drho_dt_plus_div_mass_flux_dV
relative_conservative_residual_l2
relative_last_energy_update_balance_defect
rel_operator_recon_defect
rel_output_recon_defect
```

### 7.7 MPI consistency

Run at least one relevant case with:

```bash
mpirun -np 1 ...
mpirun -np 2 ...
mpirun -np 4 ...
```

when the patch touches reductions, ownership, output, diagnostics, or parallel
I/O.

---

## 8. Current design decisions that patches must preserve

### 8.1 Option A enthalpy/species coupling

The solver uses Option A:

```text
h is the transported thermodynamic state.
T is recovered from h, Y, and p0.
```

When species transport changes `Y`, preserve transported enthalpy:

```text
h_after_species = h_before_species
T_after_species = T(h_after_species, Y_new, p0)
```

Do not reset `h` from old temperature after species changes:

```text
do not use: h = h(T_old, Y_new, p0)
```

Cantera thermo sync may combine temperature recovery and property refresh into
one pass, but it must not alter the transported `h` field.

### 8.2 Thermodynamic pressure is not projection pressure

Cantera uses:

```text
p0 = background_press
```

The flow solver pressure field is the projection pressure-like field.

Do not pass projection pressure as Cantera thermodynamic pressure unless the
governing formulation is deliberately changed and revalidated.

### 8.3 Variable-density projection target

Constant-density mode:

```text
div(u) = 0
```

Variable-density low-Mach mode:

```text
div(u) = S_projection
```

Do not use raw `div(u)` as the primary projection error in variable-density mode.

### 8.4 Conservative low-Mach source

The active source form is:

```text
S = (rho_old - rho)/(rho*dt) - (u.grad(rho))/rho
```

This targets conservative continuity:

```text
d(rho)/dt + div(rho u) = 0
```

### 8.5 Density semantics

Constant-density mode:

```text
transport%rho = params%rho
energy%rho_thermo is diagnostic
```

Variable-density mode:

```text
transport%rho <- energy%rho_thermo
```

### 8.6 Velocity BC versus mass-flow BC

Current velocity boundaries prescribe volumetric flux:

```text
given u_b -> face_flux
mass_flux = rho_face * face_flux
```

A mass-flow boundary condition is a separate future feature:

```text
given mdot -> u_n = mdot/(rho_b*A)
```

Do not silently reinterpret existing velocity boundaries as mass-flow
boundaries.

---

## 9. Documentation update rules

When a patch changes behavior, update the relevant document in the same patch or
the immediately following documentation patch.

Common mapping:

| Change type | Documents to update |
|---|---|
| New input option | `case_nml_reference.md`, `input_configuration_guide.md` |
| Numerical equation change | `numerical_method.md`, `validation_metrics.md` |
| New VTK field | `paraview_output_fields.md`, `output_layout.md` if file/layout changes |
| New CSV file/column | `output_layout.md`, `validation_metrics.md` if primary |
| Cantera bridge behavior | `cantera_interface.md`, `cantera_cpp_source_bridge.md` |
| Energy/thermo convention | `energy_thermo_enthalpy_species.md`, `numerical_method.md` |
| Build/run change | `installation_and_run.md`, `developer_guide.md` |
| FORD comment convention | `developer_guide.md` |

Do not add long patch histories to current-state user documents.  Use a changelog
or development notes for historical sequences.

---

## 10. FORD documentation rules for code patches

If a patch changes a public module, public procedure, derived type, or numerical
contract, update the FORD comments.

Use:

```text
inline math:   \( ... \)
display math:  $$ ... $$
code symbols:  `...`
```

Keep code identifiers in backticks:

```text
`transport%rho`
`energy%rho_thermo`
`fields%mass_flux`
```

Keep physical variables in math mode:

```text
\( \rho \), \( u \), \( h \), \( Y_k \), \( p_0 \)
```

Document:

```text
- units
- sign conventions
- owner-outward or cell-outward flux orientation
- time levels
- constant-density versus variable-density behavior
- whether a diagnostic is primary or explanatory
```

---

## 11. Recommended immediate next steps

The current immediate sequence should focus on hardening and clarity, not adding
chemistry yet.

Recommended order:

```text
1. Finish current-state documentation cleanup.
2. Ensure all current docs agree on:
     - active guarded variable-density mode
     - Option A enthalpy/species coupling
     - Cantera density semantics
     - projection-source versus current-source diagnostics
     - output directory layout
     - ParaView fields versus CSV validation metrics

3. Add or verify patch-wise mass-flow diagnostics for all physical boundaries.

4. Add restart read/write support for long transient-to-steady runs.

5. Centralize boundary thermodynamic state evaluation:
     h_b   = h(T_b,Y_b,p0)
     rho_b = rho(T_b,Y_b,p0)

6. Use boundary thermodynamic density to support density-aware inlet/outlet
   hardening.

7. Add optional mass-flow inlet BC only after boundary density is centralized.

8. Expand the automated validation matrix across:
     - constant-density baseline
     - Cantera-assisted constant-density
     - variable-density low-Mach
     - EOS/pressure variants
     - mesh refinement
     - dt refinement
     - MPI rank count
```

Only after the non-reacting variable-density path is robust across this matrix
should the project move to:

```text
- reaction source terms
- reaction heat release
- physical radiation coupling
- full multicomponent FV diffusion
```

---

## 12. Near-term patch candidates

### 12.1 Documentation consistency patch

Purpose:

```text
remove stale statements that variable density, species-enthalpy diffusion, or
conservative rho*h/rho*Y branches are future-only
```

Validation:

```text
FORD/docs render
manual consistency review
```

### 12.2 Boundary mass-flow diagnostics patch

Purpose:

```text
report patch-wise integrated volume flux and mass flux
```

Useful columns:

```text
patch_name
area
volume_flux_out
mass_flux_out
mean_rho_face
mean_normal_velocity
```

Validation:

```text
closed-box flow
inlet/outlet case
variable-density inlet case
MPI rank-count comparison
```

### 12.3 Restart patch

Purpose:

```text
write/read fields needed to continue a transient run
```

Minimum restart state should include:

```text
step
time
dt
u
p
face_flux or enough state to recompute it consistently
mass_flux or enough state to recompute it consistently
species Y, if enabled
enthalpy h, if enabled
temperature T, if enabled or recoverable
transport rho
transport rho_old
thermo rho_thermo
low-Mach source fields, if variable-density mode is enabled
operator-consistent rho*h state, if needed for diagnostics
```

Validation:

```text
run A continuously to time t2
run B to t1, restart, continue to t2
compare fields and diagnostics at t2
```

### 12.4 Boundary thermodynamic state helper patch

Purpose:

```text
centralize boundary h_b and rho_b evaluation from T_b,Y_b,p0
```

Validation:

```text
fixed-temperature/fixed-composition inlet
non-Dirichlet species boundary fallback to interior Y
constant-cp and Cantera thermo modes
```

### 12.5 Mass-flow boundary condition patch

Purpose:

```text
prescribe mass flux instead of velocity on selected inlet patches
```

Prerequisite:

```text
boundary density helper
```

Validation:

```text
constant-density comparison to equivalent velocity BC
variable-density inlet with known mdot
patch-wise mass-flow diagnostics
```

---

## 13. Medium-term patch candidates

```text
- bounded higher-order scalar/enthalpy convection option
- broader EOS/pressure mechanism validation cases
- HDF5/XDMF or equivalent scientific output and restart format
- improved PVTU/PVD postprocessing utilities
- manufactured qrad source tests
- physical radiation coupling through qrad
- reaction source terms and heat release
- full multicomponent diffusion interface and FV operator, if needed
```

Each of these should be split into small patches with a validation case attached
to the patch sequence.

---

## 14. Patch acceptance criteria

A patch is ready to merge when:

```text
- dry-run and apply modes work
- backups are created for changed files
- anchors are unique and robust
- git diff is narrow and intended
- debug build passes
- release build passes
- minimum validation ladder passes for affected modes
- MPI behavior is checked when relevant
- docs are updated when user-facing behavior changes
- FORD comments are updated when public code contracts change
- output schema changes are reflected in PVTU and docs
- no stale contradictory documentation remains in touched files
```

---

## 15. Maintenance rule

Keep this document focused on current patch workflow and next steps.

Do not append detailed historical patch logs here.  Historical patch sequences
belong in changelogs or development notes.

This document should answer:

```text
How should a patch be written?
How should it be applied?
How should it be validated?
Which design rules must not be broken?
What are the next development priorities?
```

