title: Installation and Run Guide


# Installation and Run Guide

This guide describes the standard setup path for building, documenting, and
running LowMachReact-Hex from a fresh checkout.

Expected workflow:

```text
1. Provide an OpenMPI compiler/runtime.
2. Create the Conda environment from environment.yml.
3. Activate the environment.
4. Build with make.
5. Run cases with mpirun through the launcher script.
6. Inspect output in ParaView and CSV diagnostics.
```

The solver is a mixed Fortran/C++ MPI application.  The Fortran solver is built
with MPI, and the C++ Cantera bridge links against the Cantera C++ library
provided by the project Conda environment.

---

## 1. Requirements

LowMachReact-Hex requires:

```text
- Conda or Mamba
- OpenMPI, including mpifort and mpirun
- a C++ compiler compatible with the Cantera C++ library
- the project Conda environment created from environment.yml
- make
```

The project Conda environment supplies:

```text
- Cantera
- FORD
- pkg-config
- Python analysis packages
- mesh and postprocessing tools
- runtime libraries used by Cantera
```

OpenMPI is normally supplied by the system or by a cluster module.  Use the same
MPI installation for building and running.

---

## 2. Load or provide OpenMPI

Before building, check that the MPI compiler and launcher are visible:

```bash
mpifort --version
mpirun --version
```

On module-based systems, load the MPI module first:

```bash
module load openmpi
```

Use the module name and version appropriate for your machine.

Important rule:

```text
The mpifort used at build time and the mpirun used at run time should come from
the same OpenMPI installation.
```

Mixing MPI installations can cause launch failures or runtime crashes.

---

## 3. Create the Conda environment

From the repository root:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate react_env
```

With Mamba:

```bash
mamba env create -f environment.yml
conda activate react_env
```

If the environment already exists, do not recreate it.  Just activate it:

```bash
conda activate react_env
```

After activation, verify:

```bash
echo "$CONDA_PREFIX"
pkg-config --modversion cantera
```

Expected:

```text
CONDA_PREFIX is non-empty
pkg-config can find Cantera
```

The Makefile checks the active Conda environment, Cantera, `pkg-config`,
`mpifort`, and runtime library paths before building.

---

## 4. Build the solver

Use the repository Makefile.

Release build:

```bash
make BUILD=release
```

Debug build with runtime checks:

```bash
make BUILD=debug
```

Convenience targets:

```bash
make release
make debug
```

The build creates:

```text
build/lowmach_react_hex.bin   compiled executable
./lowmach_react_hex           launcher script for normal runs
```

Use the launcher script for normal runs:

```bash
./lowmach_react_hex cases/channel_flow/case.nml
```

Do not run `build/lowmach_react_hex.bin` directly unless you are deliberately
debugging the build/runtime environment.  The launcher prepends the active Conda
library directory so Cantera and C++ runtime libraries are resolved from the
environment.

---

## 5. Optional build cleanup

Clean previous build products:

```bash
make clean
```

Then rebuild:

```bash
make BUILD=release
```

Use a clean rebuild after:

```text
- changing compiler or MPI modules
- recreating the Conda environment
- changing Makefile dependency order
- adding a new source module
- changing Cantera bridge build flags
```

---

## 6. `patchelf` note

If the build reports that `patchelf` is missing, install it in the active
environment:

```bash
conda install -c conda-forge patchelf
```

Then rebuild.

---

## 7. Generate meshes

Some cases include mesh-generation Makefile targets.

With the Conda environment active:

```bash
make mesh-cavity
make mesh-channel
```

These targets require the mesh tools installed by `environment.yml`, including
Gmsh-related tooling and Python mesh utilities.

The usual mesh workflow is:

```text
case geometry / mesh script
  -> native mesh files
  -> case.nml mesh_dir
  -> solver
```

Native mesh directories contain files such as:

```text
points.dat
cells.dat
faces.dat
patches.dat
periodic.dat
```

`periodic.dat` appears only for periodic connectivity.

---

## 8. Run a case

Run with OpenMPI and the launcher:

```bash
mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
```

Single-rank smoke test:

```bash
mpirun -np 1 ./lowmach_react_hex cases/lid_driven_cavity/case.nml
```

The Makefile may also provide case-specific build/run targets:

```bash
make list-cases
make channel_flow-release NP=4
make lid_driven_cavity-release NP=4
```

Use `NP` to select the number of MPI ranks for these targets.

---

## 9. Recommended smoke-test sequence

On a new machine or after changing dependencies, run a small sequence:

```bash
# 1. Debug build
make clean
make BUILD=debug

# 2. Single-rank smoke test
mpirun -np 1 ./lowmach_react_hex cases/lid_driven_cavity/case.nml

# 3. Release build
make clean
make BUILD=release

# 4. Multi-rank smoke test
mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
```

Check that:

```text
- the executable launches
- the case reads its mesh
- pressure solve does not immediately fail
- diagnostics are written when enabled
- VTK output is written when enabled
- output opens in ParaView
```

For variable-density or Cantera-heavy cases, also check that Cantera mechanism
loading and phase selection are reported correctly at startup.

---

## 10. Output files

Output location is controlled by `output_dir` in `case.nml`.

Current output layout:

```text
<output_dir>/VTK/
  ParaView VTU/PVTU/PVD files

<output_dir>/diagnostics/
  CSV diagnostics and validation files
```

Typical visualization files:

```text
*.vtu    per-rank or piece data
*.pvtu   parallel VTK master file
*.pvd    time-series collection
```

For MPI runs, open the `.pvd` or `.pvtu` file in ParaView, not a single rank
`.vtu` piece.

Typical diagnostics may include:

```text
diagnostics.csv
energy_diagnostics.csv
species_energy_conservation.csv
species_integrals.csv
enthalpy_energy_budget.csv
variable_density_diagnostics.csv
variable_density_compatibility.csv
variable_density_transport_conservation.csv
variable_density_continuity_residual.csv
```

Not every file appears in every run.  Many diagnostics are gated by enabled
physics and `write_diagnostics = .true.`.

---

## 11. ParaView quick check

After a run, open the generated `.pvd` file in ParaView to inspect the time
series.

For parallel output:

```text
Use .pvd or .pvtu for full-domain visualization.
Do not judge a multi-rank result from one rank .vtu piece.
```

Most solver fields are cell-centered finite-volume data.  ParaView line-sampling
filters may interpolate or smooth them.

---

## 12. Generate FORD documentation

FORD is installed through the Conda environment.

After activating the environment:

```bash
conda activate react_env
```

Generate documentation with the repository's FORD project file.  Typical usage:

```bash
ford ford.md
```

or, if the project file has a different name:

```bash
ford <ford-project-file>
```

FORD is primarily used for Fortran source documentation.  The C++ Cantera bridge
is documented through markdown source notes and should be linked from the docs
where practical.

Before committing documentation-related changes, check that FORD builds without
broken references or math-rendering problems.

---

## 13. Quick checklist for a new shell session

Use this when the environment already exists:

```bash
# Load MPI if your machine uses modules.
module load openmpi

# Activate project environment.
conda activate react_env

# Confirm tools.
mpifort --version
mpirun --version
pkg-config --modversion cantera

# Build or rebuild if needed.
make BUILD=release

# Run a case.
mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
```

If the environment does not exist yet, create it once:

```bash
conda env create -f environment.yml
conda activate react_env
```

---

## 14. Common setup issues

### 14.1 `CONDA_PREFIX` is empty

Activate the project environment:

```bash
conda activate react_env
```

Then verify:

```bash
echo "$CONDA_PREFIX"
```

### 14.2 `mpifort` or `mpirun` is missing

Load or install OpenMPI.

On a module system:

```bash
module load openmpi
```

Then verify:

```bash
mpifort --version
mpirun --version
```

### 14.3 Cantera is not found by `pkg-config`

Check that the project environment is active:

```bash
conda activate react_env
pkg-config --modversion cantera
```

If this still fails, recreate or repair the environment from `environment.yml`.

### 14.4 Missing `GLIBCXX` or `CXXABI` symbols at runtime

Run through the launcher with the Conda environment active:

```bash
conda activate react_env
mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
```

The launcher is designed to make the Conda C++ runtime take priority.

### 14.5 MPI mismatch

Symptoms may include launch errors, hangs, or crashes before the solver starts.

Check:

```bash
which mpifort
which mpirun
mpifort --version
mpirun --version
```

Use the same OpenMPI installation for build and run.

### 14.6 Old output confuses ParaView or scripts

If VTK arrays or output schema changed, remove the old output directory:

```bash
rm -rf cases/<case>/output
```

Then rerun the case.

### 14.7 Cantera mechanism or phase fails to load

Check the case settings:

```fortran
cantera_mech_file = "..."
cantera_phase_name = "..."
```

For a named phase, the phase name must match the YAML file exactly.

Quick Python check:

```bash
python - <<'PY'
import cantera as ct
gas = ct.Solution("mechanisms/ideal_mixavg.yaml", "gas")
print(gas.name)
print(gas.thermo_model)
print(gas.transport_model)
PY
```

---

## 15. Minimal end-to-end example

From a fresh checkout on a system with OpenMPI available:

```bash
# Load MPI if needed.
module load openmpi

# Create and activate environment.
conda env create -f environment.yml
conda activate react_env

# Verify dependencies.
mpifort --version
mpirun --version
pkg-config --modversion cantera

# Build.
make BUILD=release

# Generate meshes if needed.
make mesh-cavity
make mesh-channel

# Run.
mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
```

Then inspect:

```text
cases/channel_flow/output/VTK/
cases/channel_flow/output/diagnostics/
```

---

## 16. Maintenance notes

Keep this installation guide focused on setup and execution.

Do not add long physics explanations here.  Link to:

```text
input_configuration_guide.md
case_nml_reference.md
architecture.md
developer_guide.md
```

when users need details about runtime modes, solver physics, or validation.

Update this guide when any of the following change:

```text
- environment.yml name or dependencies
- Conda environment name
- compiler/MPI assumptions
- Makefile build targets
- launcher behavior
- mesh-generation targets
- output directory layout
- FORD documentation command
```