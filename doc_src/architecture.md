# LowMachReact-Hex Architecture

This repository implements a clean hydrodynamic foundation for reacting-flow simulations, incorporating Cantera for thermodynamics/transport and providing a scaffold for future radiation coupling.

## Mesh Pipeline

The solver utilizes a robust mesh import strategy that decouples geometry generation from the simulation core.

1.  **Gmsh Generation:** Gmsh is used to create hexahedral meshes with named physical boundary surfaces.
2.  **Native Conversion:** `tools/mesh/convert_gmsh_hex.py` reads the `.msh` file using `meshio`.
3.  **Validation:** The converter ensures all volume cells are axis-aligned cuboids (current solver requirement) and verifies connectivity.
4.  **Native Data Structures:** The converter writes native text files (`points.dat`, `cells.dat`, `faces.dat`, `patches.dat`, and `periodic.dat`).
5.  **Fortran Import:** `mod_mesh_io` reads these files into the `mesh_t` data structure.

## MPI Ownership & Decomposition

The code implements a dual-MPI architecture designed to balance flow solver efficiency with the needs of the future radiation kernel.

*   **Global Mesh Replication:** Every rank stores the full mesh (`mesh_t`). This is intentional to ensure that the radiation solver (which will decompose by wavenumber rather than space) has access to the full geometry on every rank.
*   **Cell Ownership:** `mod_mpi_flow` assigns a contiguous range of global cell IDs to each rank (`flow%first_cell` to `flow%last_cell`).
*   **Owned Updates:** Computation of cell-centered fields and momentum RHS is performed only for owned cells.
*   **Synchronization:** 
    *   **Vector Assembly:** Matrix-free operators use `flow_allreduce_global_vector` or `flow_allgather_owned_scalar` to assemble global fields from local contributions.
    *   **Scalar Reductions:** Dot products and global statistics (e.g., kinetic energy) use `MPI_Allreduce`.
*   **Separated Communicators:** 
    *   `mod_mpi_flow`: Handles flow solver and species transport synchronization.
    *   `mod_mpi_radiation`: Manages a separate communicator for radiation work items (bands, rays, external batches).

## Modular Decomposition

The solver is organized into specialized modules to isolate physics, numerics, and external dependencies:

*   **`mod_flow_projection`**: Implements the fractional-step projection algorithm.

    Advances the incompressible flow fields using the fractional-step method.

    *   Computes explicit momentum prediction (using AB2 or Forward Euler).
    *   Solves the symmetric pressure Poisson equation using matrix-free Conjugate Gradient.
    *   Corrects velocity and pressure.
    *   Calculates diagnostics (RMS divergence, net boundary flux, kinetic energy).

    > [!WARNING]
    > **Pressure Boundary Condition Limitation:** Currently, the pressure Poisson matrix is built assuming purely Neumann (zero-gradient) conditions on all physical boundaries. This gives the matrix a null space, which is removed by hardcoding `phi(1) = 0`. As a result, the solver currently only supports closed domains (cavity) or purely periodic domains (channel flow). Open flows with Dirichlet pressure inlets/outlets are structurally unsupported.

*   **`mod_species`**: Manages the transport of passive scalars ($Y_k$). Supports **Scale-on-Demand** architecture with dynamic allocation for 0 to 256+ species. Implements a **Correction Velocity** (diffusive flux correction) to ensure strict mass conservation when using different species diffusivities ($D_k$).
*   **`mod_transport_properties`**: Abstracts property evaluation. It provides a bridge to the **Cantera 3.x C++ API** for dynamic evaluation of viscosity and species diffusivity. Features include:
    *   **Automatic Mechanism Discovery**: Automatically identifies all species from a `.yaml` mechanism when `enable_reactions` is active.
    *   **Decoupled Control**: Independent toggles for Cantera-calculated fluid properties (viscosity/density) and species transport properties (diffusivity).
    *   **Name-Based Initialization**: Maps namelist-provided mass fractions to the correct indices in complex mechanisms at runtime.
*   **`mod_profiler`**: A hierarchical performance profiling module used to track execution time for critical kernels and MPI communication. It provides a terminal summary at the end of each simulation.
*   **`mod_bc`**: A unified boundary condition manager that supports field-specific types (Velocity, Pressure, Species) for every patch.

## Boundary Condition System

Boundary patches are configured in `case.nml` and mapped to mesh patches at runtime.

*   **`patch_type`**: Legacy type used as a fallback.
*   **`patch_velocity_type`**: Supports `no_slip`, `moving_wall`, `symmetry`, `periodic`.
*   **`patch_pressure_type`**: Supports `zero_gradient`, `dirichlet`, `periodic`.
*   **`patch_species_type`**: Supports `zero_gradient`, `dirichlet`, `periodic`.

## Future Growth

1.  **Energy Equation:** Coupling temperature/enthalpy into the projection method.
2.  **Variable Density:** Transitioning from incompressible to the low-Mach formulation.
3.  **Chemistry Source Terms:** Integrating Cantera reaction rates into `mod_species`.
4.  **Radiation Integration:** Coupling the wavenumber-decomposed `q_rad` source term into the energy equation.
