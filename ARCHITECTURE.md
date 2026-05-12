# Architecture Overview (LowMachReact-Hex)

This document describes the high-level system design, module responsibilities, and data flow of the LowMachReact-Hex solver.

## System Design Philosophy

The solver is designed for **Research and Scalability**. It prioritizes physical accuracy and developer transparency over "black-box" performance.

### Key Architectural Principles:
1.  **Replicated Mesh**: Every MPI rank holds the full topological structure of the mesh. This is a non-negotiable requirement to support future spectral radiation solvers that need global geometric visibility.
2.  **Owned-Cell Decomposition**: While the mesh is replicated, the **computational work** for the flow and species fields is partitioned. Each rank is responsible for a contiguous range of "owned" cells.
3.  **Scale-on-Demand Species**: The species transport system is dynamic. It can handle 0 to 256+ species without recompilation, using automatic discovery from Cantera mechanism files.
4.  **Decoupled Physics**: The flow solver, species transport, and property evaluation (Cantera) are decoupled. This allows for easy testing of individual components.

---

## Module Responsibilities

### Core Modules
| Module | Responsibility |
| :--- | :--- |
| `mod_kinds` | Precision definitions and global constants. |
| `mod_input` | Namelist parsing and case configuration. |
| `mod_mesh_types` | Geometric data structures (cells, faces, patches). |
| `mod_mesh_io` | Native binary mesh ingestion. |
| `mod_fields` | Storage for velocity, pressure, and species fields. |

### Physics & Numerics
| Module | Responsibility |
| :--- | :--- |
| `mod_flow_projection` | Fractional-step incompressible Navier-Stokes solver. |
| `mod_species` | Advection-diffusion transport for mass fractions. |
| `mod_bc` | Boundary condition evaluation and aliasing. |
| `mod_transport_properties` | Interface to Cantera for $\mu, \rho, D_k$. |

### Parallelism & Infrastructure
| Module | Responsibility |
| :--- | :--- |
| `mod_mpi_flow` | Domain decomposition and collective communication. |
| `mod_profiler` | Hierarchical performance monitoring. |
| `mod_output` | VTU/PVD and diagnostic CSV generation. |

---

## Data Flow

### Timestep Execution Path:
1.  **Predictor**: `mod_flow_projection` computes an intermediate velocity $\mathbf{u}^*$ using AB2 or Forward Euler.
2.  **Poisson RHS**: `mod_flow_projection` computes the divergence of the predicted face fluxes.
3.  **Poisson Solve**: PCG solver (in `mod_flow_projection`) finds the pressure potential $\phi$.
4.  **Correction**: Velocity and pressure fields are updated to satisfy the divergence-free constraint.
5.  **Species Transport**: `mod_species` advances the mass fractions using the corrected fluxes.
6.  **Diagnostics**: `mod_output` and `mod_profiler` generate status reports.

---

## MPI Communication Model

LowMachReact-Hex uses a hybrid communication model:
-   **Halo Exchanges**: Used for neighbor-cell access during flux and gradient calculations.
-   **Global Gathers**: `MPI_Allgatherv` is used to synchronize owned-cell updates across the replicated mesh.
-   **Reductions**: `MPI_Allreduce` is used for global statistics (kinetic energy, max divergence) and PCG dot products.
