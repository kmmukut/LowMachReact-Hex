# LowMachReact-Hex Solver

This repository contains a Fortran 2008 MPI finite-volume solver designed for laminar incompressible and low-Mach-style flows on hexahedral meshes. The solver is built for high-performance physics simulations, with a focus on extensions into reacting flows and radiation heat transfer.

## Project Overview

The current baseline implementation features:
* **Cell-centered finite-volume** formulation on axis-aligned cuboid cells.
* **Fractional-step projection method** with corrected conservative face fluxes.
* **MPI decomposition** by owned cell ranges, with a globally replicated mesh to support spectral radiation solvers.
* **Scale-on-Demand Species**: Dynamic passive species transport from 0 to 256+ species. Implements diffusive-flux correction for strict mass-fraction consistency.
* **Passive Sensible-Enthalpy Energy Transport**: Optional transport of sensible enthalpy `h`, with temperature recovered from `h`, composition `Y`, and thermodynamic pressure `p0`.
* **Cantera 3.x Integration**: Decoupled Cantera paths for flow transport, species diffusivity, and energy thermodynamics. Cantera thermo provides `h <-> T`, `cp`, thermal conductivity, and `rho_thermo`.
* **Guarded Variable-Density Mode**: Experimental non-reacting low-Mach mode can sync active density from Cantera thermo density for the supported Cantera-energy configuration.
* **Integrated Profiling**: Hierarchical MPI-aware profiler for monitoring projection, species, energy, Cantera thermo sync, I/O, and communication overhead.
* **VTU/PVD output** for visualization in ParaView and CSV diagnostics for quantitative monitoring.

## Governing Equations

The stable baseline is the incompressible, constant-density, laminar
Navier-Stokes system:

$$
\nabla \cdot \mathbf{u} = 0
$$

$$
\frac{\partial \mathbf{u}}{\partial t} + \nabla \cdot (\mathbf{u} \mathbf{u}) = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{f}_{body}
$$

The guarded variable-density low-Mach path changes the projection target to:

$$
\nabla \cdot \mathbf{u} = S_{\mathrm{projection}}
$$

with active density from the selected Cantera phase when
`enable_variable_density=.true.` and `density_eos="cantera"`.

Passive, non-reacting species transport is also supported. In constant-density
mode it transports \(Y_k\); in variable-density mode it uses conservative
\(\rho Y_k\):

$$
\frac{\partial Y_k}{\partial t} + \nabla \cdot (\mathbf{u} Y_k) = \nabla \cdot (D_k \nabla Y_k)
$$

Optional sensible-enthalpy transport is available:

$$
\frac{\partial h}{\partial t}
+
\nabla \cdot (\mathbf{u}h)
=
\frac{1}{\rho}\nabla \cdot (\lambda \nabla T)
+
\frac{q_{rad}}{\rho}
$$

The transported thermal state is `h`; temperature is recovered from current
`h`, current `Y`, and uniform thermodynamic pressure `p0`. In constant-density
mode `transport%rho = params%rho` and Cantera `rho_thermo` is diagnostic. In
guarded variable-density mode `transport%rho <- energy%rho_thermo`, so the
actual EOS is the selected Cantera YAML phase, not the string `"cantera"`.

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
*   **[Developer Guide](doc_src/DEVELOPER_GUIDE.md)**: Development workflow, FORD conventions, and patch discipline.
*   **[Architecture](doc_src/architecture.md)**: Detailed split of responsibilities and data structures.
*   **[Numerical Method](doc_src/numerics.md)**: Equations, discretization, time ordering, and implemented physics.
*   **[Input Reference](doc_src/input_configuration.md)** and **[Input Guide](doc_src/input_configuration_guide.md)**: Compact namelist reference and detailed configuration guidance.
*   **[Profiling Guide](doc_src/profiling.md)**: How to read the nested MPI-aware profiling report.

To generate the full API documentation using FORD:
```bash
ford  ford.md
```
