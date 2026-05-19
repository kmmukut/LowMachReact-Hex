# LowMachReact-Hex Documentation

Welcome to the technical documentation for **LowMachReact-Hex**, a powerful agentic CFD solver for low-Mach number reacting flows on hexahedral meshes.

## Documentation Sections

- [Numerical Method](numerics.md): Detailed description of the governing equations, time integration, and spatial discretization.
- [Architecture](architecture.md): Overview of the system design, MPI decomposition, and module structure.
- [Input Configuration](input_configuration.md): Guide to the `case.nml` parameters and boundary condition setup.

## Project Goal

The primary objective of this solver is to provide a scalable, research-oriented framework for simulating complex multi-species transport with future extensions for spectral radiation and finite-rate chemistry.

### Key Features:
- **Fractional-Step Projection**: Robust pressure-velocity coupling for incompressible flows.
- **Scale-on-Demand Species**: Dynamic handling of arbitrary species counts.
- **Replicated Mesh MPI**: Optimized for global geometric visibility, required for future radiation solvers.
- **Cantera Integration**: State-of-the-art thermodynamics and transport property evaluation.

---
*Generated using FORD (Fortran Documenter).*
