# Architecture Overview (LowMachReact-Hex)

This document describes the high-level system design, module responsibilities, and data flow of the LowMachReact-Hex solver.

## System Design Philosophy

The solver is designed for **Research and Scalability**. It prioritizes physical accuracy and developer transparency over "black-box" performance.

### Key Architectural Principles:
1.  **Replicated Mesh**: Every MPI rank holds the full topological structure of the mesh. This is a non-negotiable requirement to support future spectral radiation solvers that need global geometric visibility.
2.  **Owned-Cell Decomposition**: While the mesh is replicated, the **computational work** for the flow and species fields is partitioned. Each rank is responsible for a contiguous range of "owned" cells.
3.  **Scale-on-Demand Species**: The species transport system is dynamic. It can handle 0 to 256+ species without recompilation, using explicit species lists and Cantera species-name mapping when enabled.
4.  **Passive Enthalpy as Energy State**: The energy path transports sensible enthalpy `h`; temperature and thermodynamic properties are recovered from `h`, `Y`, and `p0`.
5.  **Decoupled Physics**: The flow solver, species transport, energy transport, and property evaluation are decoupled. This allows for focused validation of individual components.

---

## Module Responsibilities

### Core Modules
| Module | Responsibility |
| :--- | :--- |
| `mod_kinds` | Precision definitions and global constants. |
| `mod_input` | Namelist parsing and case configuration. |
| `mod_mesh_types` | Geometric data structures (cells, faces, patches). |
| `mod_mesh_io` | Native text mesh ingestion. |
| `mod_fields` | Storage for velocity, pressure, and conservative face-flux fields. |
| `mod_energy` | Storage and update routines for temperature, sensible enthalpy, `qrad`, `cp`, thermal conductivity, and Cantera `rho_thermo`. |

### Physics & Numerics
| Module | Responsibility |
| :--- | :--- |
| `mod_flow_projection` | Fractional-step projection solver: constant-density \( \nabla\cdot u=0 \) baseline plus guarded variable-density \( \nabla\cdot u=S \). |
| `mod_species` | Non-reacting species transport for `Y` or conservative `rho*Y` in variable-density mode. |
| `mod_bc` | Boundary condition evaluation and aliasing. |
| `mod_transport_properties` | Interface to constant or Cantera transport properties such as \(\mu\) and \(D_k\); owns active `transport%rho`. |

### Parallelism & Infrastructure
| Module | Responsibility |
| :--- | :--- |
| `mod_mpi_flow` | Domain decomposition and collective communication. |
| `mod_profiler` | Hierarchical MPI-aware performance monitoring, including projection, species, energy, Cantera thermo sync, output, and MPI communication regions. |
| `mod_output` | VTU/PVD and diagnostic CSV generation. |

---

## Data Flow

### Timestep Execution Path:
1.  **Transport-property refresh**: `mod_transport_properties` updates Cantera or constant transport properties on `transport_update_interval`.
2.  **Predictor**: `mod_flow_projection` computes an intermediate velocity $\mathbf{u}^*$ using AB2 or Forward Euler.
3.  **Poisson RHS**: `mod_flow_projection` computes the divergence of the predicted face fluxes.
4.  **Poisson Solve**: PCG solver in `mod_flow_projection` finds the pressure potential $\phi$.
5.  **Correction**: Velocity, projection pressure, volumetric face flux, and density-weighted face flux are updated.
6.  **Species Transport**: `mod_species` advances non-reacting mass fractions using volumetric flux in constant-density mode or mass flux for conservative `rho*Y` in variable-density mode.
7.  **Energy Transport**: `mod_energy` advances sensible enthalpy. With Cantera thermo, dependent state is synchronized as `(T, cp, lambda, rho_thermo) = sync(h,Y,p0)`; variable-density mode then syncs `transport%rho <- energy%rho_thermo`.
8.  **Diagnostics and Output**: `mod_output` and `mod_profiler` generate CSV, VTU/PVTU/PVD, and terminal status reports.

---

## MPI Communication Model

LowMachReact-Hex uses a hybrid communication model:
-   **Halo Exchanges**: Used for neighbor-cell access during flux and gradient calculations.
-   **Global Gathers**: `MPI_Allgatherv` is used to synchronize owned-cell updates across the replicated mesh.
-   **Reductions**: `MPI_Allreduce` is used for global statistics, validation metrics, and PCG dot products. In variable-density mode, projection validation uses `div(u)-S_projection` rather than raw divergence.
