# LowMachReact-Hex Solver

This repository contains a Fortran 2008 MPI finite-volume solver designed for laminar incompressible and low-Mach-style flows on hexahedral meshes. The solver is built for high-performance physics simulations, with a focus on extensions into reacting flows and radiation heat transfer.

## Project Overview

The current baseline implementation features:
* **Cell-centered finite-volume** formulation on axis-aligned cuboid cells.
* **Fractional-step projection method** with corrected conservative face fluxes.
* **MPI decomposition** by owned cell ranges, with a globally replicated mesh to support spectral radiation solvers.
* **Scale-on-Demand Species**: Dynamic species transport from 0 to 256+ species. Implements **Correction Velocity** for strict mass conservation. Supports manual selection or automatic discovery from Cantera mechanisms.
* **Cantera 3.x integration**: Decoupled flags for viscosity/density and species diffusivity. Automatic mechanism loading for reacting flows.
* **Integrated Profiling**: Hierarchical MPI-aware profiler for monitoring kernel execution time and communication overhead.
* **VTU/PVD output** for visualization in ParaView and CSV diagnostics for quantitative monitoring.

## Governing Equations

The solver represents incompressible, constant-density, laminar Navier-Stokes equations:

$$
\nabla \cdot \mathbf{u} = 0
$$

$$
\frac{\partial \mathbf{u}}{\partial t} + \nabla \cdot (\mathbf{u} \mathbf{u}) = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{f}_{body}
$$

Passive species transport is also supported:

$$
\frac{\partial Y_k}{\partial t} + \nabla \cdot (\mathbf{u} Y_k) = \nabla \cdot (D_k \nabla Y_k)
$$

## Installation & Environment

The project requires a Fortran 2008 compiler and the Cantera 3.x C++ library. It is recommended to use the provided Conda environment:

```bash
# 1. Create the environment from the provided YAML file
conda env create -f environment.yml

# 2. Activate the environment
conda activate react_env
```

## Quick Start

Once the environment is active:

```bash
# 1. Build the solver
make BUILD=release

# 2. Run a validation case (e.g., channel flow)
mpirun -np 8 ./lowmach_react_hex cases/channel_flow/case.nml
```

Visualization files are written to the `output/` directory of the case. Open `flow.pvd` in ParaView to view the time series.

## Documentation

*   **[Live Documentation](https://kmmukut.github.io/LowMachReact-Hex)**: Full API reference, call graphs, and search.
*   **[Developer Guide](DEVELOPER_GUIDE.md)**: Strict rules, MPI architecture, and non-negotiable design principles.
*   **[Architecture](doc_src/architecture.md)**: Detailed split of responsibilities and data structures.
*   **[Numerics](doc_src/numerics.md)**: Mathematical discretization and solver details.
*   **[Input Configuration](doc_src/input_configuration.md)**: Detailed guide on namelist parameters and input file preparation.

To generate the full API documentation using FORD:
```bash
ford ford.yml
```
