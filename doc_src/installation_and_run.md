title: Installation and Run Guide

# Installation and Run Guide

This page describes the standard setup path for building and running
LowMachReact-Hex from a fresh checkout. The expected workflow is:

1. Provide an OpenMPI compiler/runtime.
2. Create the Conda environment from `environment.yml`.
3. Activate the environment.
4. Build with `make`.
5. Run a case with `mpirun`.

## Requirements

LowMachReact-Hex requires:

- Conda or Mamba for the project environment.
- OpenMPI, including `mpifort` and `mpirun`.
- A C++ compiler compatible with the Cantera C++ library.
- The project Conda environment created from `environment.yml`.

The Conda environment supplies Cantera, FORD, mesh tooling, Python analysis
packages, and `pkg-config`. OpenMPI is normally supplied by the system or by a
cluster module. Use the same MPI installation for building and running.

Check that OpenMPI is visible before building:

```bash
mpifort --version
mpirun --version
```

On module-based systems this may require a command such as:

```bash
module load openmpi
```

Use the module name/version appropriate for your machine.

## Create the Conda Environment

From the repository root, create the project environment:

```bash
conda env create -f environment.yml
```

Then activate it:

```bash
conda activate react_env
```

If you prefer Mamba, the equivalent command is:

```bash
mamba env create -f environment.yml
conda activate react_env
```

After activation, confirm that the environment is active and that Cantera is
visible through `pkg-config`:

```bash
echo "$CONDA_PREFIX"
pkg-config --modversion cantera
```

The Makefile checks `CONDA_PREFIX`, Cantera, `pkg-config`, `mpifort`, and the
runtime library path before building.

## Build the Solver

For a release build:

```bash
make BUILD=release
```

For a debug build with runtime checks:

```bash
make BUILD=debug
```

Convenience targets are also available:

```bash
make release
make debug
```

The build creates:

- `build/lowmach_react_hex.bin`: the compiled executable.
- `./lowmach_react_hex`: a launcher script used for normal runs.

Run the launcher, not the binary directly. The launcher prepends the active
Conda library directory so the Cantera and C++ runtime libraries come from the
environment.

If the build reports that `patchelf` is missing, install it in the active
environment:

```bash
conda install -c conda-forge patchelf
```

## Generate Meshes

Some cases include mesh-generation scripts. With the Conda environment active:

```bash
make mesh-cavity
make mesh-channel
```

These targets require the mesh tools from `environment.yml`, including `gmsh`
and `meshio`.

## Run a Case

Run the solver with OpenMPI:

```bash
mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
```

For a single-rank smoke test:

```bash
mpirun -np 1 ./lowmach_react_hex cases/lid_driven_cavity/case.nml
```

The Makefile can also build and run known cases. Use `NP` to choose the number
of MPI ranks:

```bash
make list-cases
make channel_flow-release NP=4
make lid_driven_cavity-release NP=4
```

## Output Files

Simulation output is written according to the `output_dir` value in the case
`case.nml`. Typical outputs include:

- VTU/PVTU files for ParaView.
- PVD collection files for time-series loading.
- CSV diagnostics for flow, energy, and profiling-related analysis.

Open the generated `.pvd` file in ParaView to inspect the time series.

## Quick Checklist

Use this checklist when setting up a new machine or shell session:

```bash
# Load or otherwise provide OpenMPI.
module load openmpi

# Create once, then activate each session.
conda env create -f environment.yml
conda activate react_env

# Build.
make BUILD=release

# Run.
mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
```

If the Conda environment already exists, skip the `conda env create` step and
only run `conda activate react_env`.

## Common Setup Issues

If `CONDA_PREFIX` is empty, activate the environment:

```bash
conda activate react_env
```

If `mpifort` or `mpirun` is missing, load or install OpenMPI.

If `pkg-config --modversion cantera` fails, confirm that `react_env` was created
from the repository `environment.yml` and is currently active.

If the executable starts but reports missing `GLIBCXX` or `CXXABI` symbols, run
through `./lowmach_react_hex` with the Conda environment active. The launcher is
designed to make the Conda C++ runtime take priority.
