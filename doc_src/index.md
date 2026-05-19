title: LowMachReact-Hex Documentation

# LowMachReact-Hex Documentation

Welcome to the technical documentation for **LowMachReact-Hex**, a Fortran 2008 MPI finite-volume solver for projection/low-Mach flow on hexahedral meshes.

The stable baseline is constant-density incompressible projection. The current guarded development path supports non-reacting variable-density low-Mach coupling with Cantera thermo density, conservative `rho*Y` species transport, conservative `rho*h` enthalpy transport, VTU/PVTU output, validation diagnostics, and MPI-aware profiling. Reactions, reaction heat release, mass-flow boundary conditions, restart, and external radiation physics remain future work.

## Documentation Sections

- [Numerical Method](numerics.html): Governing equations, finite-volume discretization, projection method, species transport, and sensible-enthalpy transport.
- [Current Solver State](current_solver_state.html): Current enabled physics, disabled physics, runtime ordering, and density/thermo conventions.
- [Installation and Run Guide](installation_and_run.html): Conda environment creation, OpenMPI setup, build commands, mesh generation, and case execution.
- [Input Configuration](input_configuration.html): Compact `case.nml` reference.
- [Input Configuration Guide](input_configuration_guide.html): Longer configuration guide with recommended validation modes and examples.
- [Energy and Thermodynamic Conventions](energy_thermo_conventions.html): Sensible enthalpy, temperature recovery, Cantera reference state, heat conduction, and `qrad` sign convention.
- [Cantera Interface Notes](cantera_interface.html): Cantera bridge responsibilities, sensible enthalpy convention, composition rules, density rules, and cache dependencies.
- [Output Layout](output_layout.html): Solver-generated directories and files.
- [ParaView Output Fields](paraview_output_fields.html): Field meanings and ParaView interpretation.
- [Validation Metrics](validation_metrics.html): Primary pass/fail metrics and interpretation.
- [Validation Automation](validation_automation.html): Checker, audit, and matrix-runner workflow.
- [Profiling Guide](profiling.html): How to enable and interpret the MPI-aware nested profiler.
- [Architecture](architecture.html): System design, MPI decomposition, module responsibilities, and future growth path.

- [DEVELOPER_GUIDE.](DEVELOPER_GUIDE.html): Instructions for development.

## Current Solver Stage

The current solver has a stable constant-density baseline:

```text
rho_flow = params%rho
p0       = params%background_press
h        = transported sensible enthalpy
T        = T(h, Y, p0)
```

In this mode, Cantera thermodynamic density is available as `rho_thermo` when
thermo is enabled, but it is diagnostic and is not used by the pressure
projection.

The guarded variable-density path uses:

```text
transport%rho <- energy%rho_thermo
projection target = div(u) = S_projection
species/energy = conservative rho*Y and rho*h
```

For that path, raw `div(u)` is explanatory; projection validation uses
`divu_minus_S_projection_*`.

## Project Goal

The primary objective is to provide a scalable, research-oriented framework for validating finite-volume flow, species, enthalpy, thermodynamic-property, and future radiation/chemistry coupling on hexahedral meshes.

## Key Features

- **Fractional-Step Projection**: Pressure-velocity coupling for constant-density incompressible flow.
- **Scale-on-Demand Species**: Dynamic handling of passive species mass fractions.
- **Passive Sensible-Enthalpy Transport**: Transported `h` with temperature recovered from `h`, `Y`, and `p0`.
- **Cantera Integration**: Transport properties and thermodynamics through Cantera, including `h <-> T`, `cp`, `lambda`, and `rho_thermo`.
- **Replicated Mesh MPI**: Full mesh visibility on each rank with owned-cell decomposition for current flow/species/energy work.
- **MPI-Aware Profiling**: Inclusive flat and nested timing reports for projection, species, energy, Cantera sync, output, and MPI communication.
- **VTU/PVTU Output**: ParaView-ready visualization output and CSV diagnostics.

---
*Generated using FORD (Fortran Documenter).*
