title: LowMachReact-Hex Documentation

# LowMachReact-Hex Documentation

Welcome to the technical documentation for **LowMachReact-Hex**, a Fortran 2008 MPI finite-volume solver for constant-density low-Mach / incompressible flow on hexahedral meshes.

The current staged solver supports incompressible projection, passive species transport, passive sensible-enthalpy transport, Cantera thermodynamic/property evaluation, VTU/PVTU output, and MPI-aware profiling. Variable-density low-Mach coupling, reactions, reaction heat release, and external radiation physics remain future work.

## Documentation Sections

- [Numerical Method](numerics.html): Governing equations, finite-volume discretization, projection method, species transport, and sensible-enthalpy transport.
- [Current Solver State](current_solver_state.html): Current enabled physics, disabled physics, runtime ordering, and density/thermo conventions.
- [Installation and Run Guide](installation_and_run.html): Conda environment creation, OpenMPI setup, build commands, mesh generation, and case execution.
- [Input Configuration](input_configuration.html): Compact `case.nml` reference.
- [Input Configuration Guide](input_configuration_guide.html): Longer configuration guide with recommended validation modes and examples.
- [Energy and Thermodynamic Conventions](energy_thermo_conventions.html): Sensible enthalpy, temperature recovery, Cantera reference state, heat conduction, and `qrad` sign convention.
- [Enthalpy/Species Coupling Convention](enthalpy_species_coupling_convention.html): Option A convention for preserving transported `h` when species composition changes.
- [Cantera Interface Notes](cantera_interface.html): Cantera bridge responsibilities, sensible enthalpy convention, composition rules, density rules, and cache dependencies.
- [Profiling Guide](profiling.html): How to enable and interpret the MPI-aware nested profiler.
- [Validation Ladder](validation_ladder.html): Recommended build, regression, energy, Cantera thermo, species, and future `qrad` validation steps.
- [Architecture](architecture.html): System design, MPI decomposition, module responsibilities, and future growth path.

- [DEVELOPER_GUIDE.](DEVELOPER_GUIDE.html): Instructions for development.

## Current Solver Stage

The current solver should be interpreted as a constant-density projection code with optional passive scalar physics:

```text
rho_flow = params%rho
p0       = params%background_press
h        = transported sensible enthalpy
T        = T(h, Y, p0)
```

Cantera thermodynamic density is available as `rho_thermo`, but it is diagnostic/future-use only and is not used by the pressure projection.

## Project Goal

The primary objective is to provide a scalable, research-oriented framework for validating finite-volume flow, species, enthalpy, thermodynamic-property, and future radiation/chemistry coupling on hexahedral meshes.

## Key Features

- **Fractional-Step Projection**: Pressure-velocity coupling for constant-density incompressible flow.
- **Scale-on-Demand Species**: Dynamic handling of passive species mass fractions.
- **Passive Sensible-Enthalpy Transport**: Transported `h` with temperature recovered from `h`, `Y`, and `p0`.
- **Cantera Integration**: Transport properties and thermodynamics through Cantera, including `h <-> T`, `cp`, `lambda`, and diagnostic `rho_thermo`.
- **Replicated Mesh MPI**: Full mesh visibility on each rank with owned-cell decomposition for current flow/species/energy work.
- **MPI-Aware Profiling**: Inclusive flat and nested timing reports for projection, species, energy, Cantera sync, output, and MPI communication.
- **VTU/PVTU Output**: ParaView-ready visualization output and CSV diagnostics.

---
*Generated using FORD (Fortran Documenter).*
