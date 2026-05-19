# Non-reacting counterflow (Gmsh)

Structured hexahedral brick with **opposed axial inlets** on `xmin` and `xmax`.
The mesh is a **thin 2D slab**: resolved in `x-y`, exactly **one hex cell
thick in `z`** (`nz = 1` in `counterflow.geo`). The transverse `ymin` / `ymax`
patches are open `neumann` outlets so the constant-density projection solve can
discharge the mass injected by the opposed inlets. The `zmin` / `zmax` patches
are `symmetric`, giving a one-cell-thick slab rather than a full 3D duct.

Patch names follow `cases/channel_flow/channel.geo` so they map directly to `case.nml`.

## Environment

From the repository root, load modules then activate the conda env ( **`gmsh` is on `PATH` from `react_env`** after activation):

```bash
module load openmpi miniforge3 git
conda activate react_env
```

Headless automation (e.g. default sandbox agents) typically **cannot** run that activation sequence, so `gmsh`, `meshio`, and the solver binary are **not** assumed to exist there. Run meshing, the executable, and the Python tools **in your own terminal** with `react_env` active, as you did.

## Mesh

```bash
cd cases/counterflow_nonreacting_gmsh
./make_mesh.sh
```

This runs `gmsh -3 counterflow.geo` and `tools/mesh/convert_gmsh_hex.py` to write `mesh_native/`. No periodic pairs are passed to the converter (unlike the channel case).

The solver uses **3D hexes only**; this case therefore uses a **degenerate-thickness 3D slab** (`nz = 1`) instead of a native 2D grid. If you previously meshed the old quasi-1D geometry, rerun `./make_mesh.sh` so the new `64 x 32 x 1` mesh is generated.

## Run

From the project root, after building the solver:

```bash
mpirun -np 1 ./lowmach_react_hex cases/counterflow_nonreacting_gmsh/case.nml
```

## Boundary intent

| Patch   | Role            | Velocity (`u,v,w`)   | Species / `T`        |
|---------|-----------------|------------------------|------------------------|
| `xmin`  | Fuel-side inlet | `(+U0, 0, 0)`          | Dirichlet (see nml)  |
| `xmax`  | Oxidizer inlet  | `(-U0, 0, 0)`          | Dirichlet            |
| `ymin`  | Open outlet     | zero-gradient / balanced flux | Zero-gradient scalars  |
| `ymax`  | Open outlet     | zero-gradient / balanced flux | Zero-gradient scalars  |
| `zmin`  | Symmetry plane  | slip symmetry          | Zero-gradient scalars  |
| `zmax`  | Symmetry plane  | slip symmetry          | Zero-gradient scalars  |

Default `U0 = 0.5` m/s in `case.nml`; adjust with the Cantera 1D reference strain / inlet velocity for validation.

## Post-run validation (FV centerline + optional Cantera 1D)

With `react_env` active (`numpy`, `pandas`, `meshio`, `matplotlib`, `cantera` as in `environment.yml`), from the repo root.

### 1) Cantera non-reacting 1D reference (same inlets / width as this case)

GRI-Mech 3.0 is **bundled with Cantera** as ``gri30.yaml``; the script default loads it from Cantera’s data directory (no copy in this repo required). Use ``--mech /path/to/foo.yaml`` only for a custom mechanism.

```bash
python tools/cantera_counterflow_reference.py \
  --output cases/counterflow_nonreacting_gmsh/validation/cantera_reference.csv
```

This builds a `CounterflowDiffusionFlame`, sets **all reaction multipliers to zero**, matches **fuel/oxidizer mass fractions** and **Lx = 2 m** to `case.nml` / `counterflow.geo`, and writes CSV columns `x`, `u`, `temperature`, `Y_CH4`, `Y_O2`, `Y_N2`.

### 2) Extract the FV centerline from VTU (merge MPI pieces via `.pvtu`)

```bash
python tools/validate_flow.py counterflow \
  --vtu cases/counterflow_nonreacting_gmsh/output/flow_050000.pvtu
```

Writes `validation/counterflow_centerline.csv`, `counterflow_summary.txt`, and **`counterflow_centerline.png`** by default (`--no-plot` skips the PNG). Matplotlib uses a non-interactive backend so the figure saves without a display.

### 3) Overlay errors vs the Cantera CSV (interpolates FV onto reference `x`)

```bash
python tools/validate_flow.py counterflow \
  --vtu cases/counterflow_nonreacting_gmsh/output/flow_050000.pvtu \
  --reference-csv cases/counterflow_nonreacting_gmsh/validation/cantera_reference.csv
```

This writes `validation/counterflow_vs_reference_errors.csv` for every column name present in **both** the reference and the FV profile (e.g. `u`, `temperature`, `Y_CH4`). With plotting on (default), it also updates **`validation/counterflow_centerline.png`** (FV vs Cantera overlaid) and writes **`validation/counterflow_vs_reference_residuals.png`** (FV−Cantera on the reference `x` grid, up to four quantities).

## Validation note

Treat this case as a **counterflow validation harness**, not yet as a passing
quantitative benchmark. Cantera's `CounterflowDiffusionFlame` is a 1D similarity
solution with opposed-flow stagnation behavior, transverse/radial strain, and
variable thermodynamic density. The current FV solver is still documented as a
constant-density projection solver with conservative face fluxes based on the
configured `rho` and `nu`.

For the present code, use this case to check:

- the Cantera reference can be regenerated reproducibly,
- `Y_k` remains bounded and `sum_Y` remains close to one,
- the FV extraction/overlay tooling works,
- diagnostics expose whether the flow field is actually compatible with a
  counterflow comparison.

A credible Cantera comparison should show a centerline axial velocity that
changes sign between the two inlets and small divergence/mass-balance errors.
If `validate_flow.py counterflow --reference-csv ...` reports no `u` sign change
or large divergence, do **not** interpret the FV-vs-Cantera error norms as
validation pass/fail numbers. They are telling you the flow model/configuration
is not yet the same problem as Cantera's opposed-flow test.
