# Input Configuration Guide

This document provides a detailed explanation of how to configure simulations for the LowMachReact-Hex solver using the `case.nml` namelist and preparing the necessary input files.

## 1. The `case.nml` File

The solver uses a standard Fortran Namelist format. The file is divided into several blocks:

### `&mesh_input`
*   **`mesh_dir`**: Path to the directory containing the native mesh files (`points.dat`, `cells.dat`, `faces.dat`, etc.).

### `&time_input`
*   **`nsteps`**: Total number of timesteps to run.
*   **`dt`**: Initial timestep size (seconds).
*   **`use_dynamic_dt`**: Boolean (`.true.`/`.false.`). If true, the solver adjusts `dt` based on the target `max_cfl`.
*   **`max_cfl`**: The target maximum Courant-Friedrichs-Lewy number for dynamic time-stepping.
*   **`output_interval`**: How often (in steps) to write VTU/PVD files and diagnostic rows.

### `&fluid_input`
*   **`rho`**: Constant density ($kg/m^3$).
*   **`nu`**: Constant kinematic viscosity ($m^2/s$). Only used if `enable_cantera` is false.
*   **`enable_cantera`**: Boolean. If true, fluid viscosity is calculated dynamically using a Cantera mechanism.
*   **`cantera_mech_file`**: Path to the Cantera `.yaml` or `.xml` mechanism file.
*   **`transport_update_interval`**: How often (in steps) to recalculate transport properties (viscosity, diffusivity).
*   **`background_temp`**: Temperature used for property calculation if the energy equation is absent.
*   **`background_press`**: Pressure used for property calculation.

### `&solver_input`
*   **`pressure_max_iter`**: Maximum Conjugate Gradient iterations for the pressure Poisson solver.
*   **`pressure_tol`**: Convergence tolerance for the pressure solver.
*   **`convection_scheme`**: Either `"upwind"` (1st order, stable) or `"central"` (2nd order, accurate but less stable).
*   **`body_force_x/y/z`**: Volumetric body forces (e.g., gravity).

### `&boundary_input`
*   **`n_patches`**: Number of boundary patches defined in the mesh.
*   **`patch_name`**: List of strings matching the names in `patches.dat`.
*   **`patch_velocity_type`**: `"no_slip"`, `"slip"`, `"dirichlet"`, `"zero_gradient"`, or `"periodic"`.
*   **`patch_pressure_type`**: `"dirichlet"`, `"zero_gradient"`, or `"periodic"`.
*   **`patch_species_type`**: `"dirichlet"`, `"zero_gradient"`, or `"periodic"`.
*   **`patch_u/v/w/p`**: Specified values for Dirichlet boundaries.
*   **`patch_Y(species_id, patch_id)`**: Specified mass fractions for species boundaries.

### `&species_input`
*   **`enable_species`**: Master toggle for species transport.
*   **`nspecies`**: Number of species to transport.
*   **`species_name`**: List of species names (e.g., `"O2"`, `"N2"`).
*   **`species_diffusivity`**: Constant diffusivity coefficients ($m^2/s$) if Cantera is disabled.
*   **`initial_Y`**: Initial mass fraction distribution in the domain.
*   **`enable_cantera`**: Boolean. If true, species diffusivities are calculated dynamically via Cantera.

### `&output_input`
*   **`output_dir`**: Directory where VTU and CSV files will be saved.
*   **`write_vtu`**: Boolean. Enable/disable volume output.
*   **`write_diagnostics`**: Boolean. Enable/disable `diagnostics.csv`.

---

## 2. Preparing Input Files

### Native Mesh Format
The solver expects a directory containing the following ASCII files:

1.  **`points.dat`**:
    *   Header: `npoints`
    *   Lines: `id x y z`
2.  **`cells.dat`**:
    *   Header: `ncells`
    *   Lines: `id node1 node2 ... node8 cx cy cz volume`
3.  **`faces.dat`**:
    *   Header: `nfaces`
    *   Lines: `id owner neighbor patch nx ny nz area cx cy cz`
4.  **`patches.dat`**:
    *   Header: `npatches`
    *   For each patch: `id name nfaces` followed by a list of face IDs.
5.  **`periodic.dat`** (Optional):
    *   Header: `nlinks`
    *   Lines: `face_id pair_face_id neighbor_cell_id`

### Cantera Mechanisms
If Cantera is enabled, you must provide a valid mechanism file (typically `gri30.yaml`). The solver will automatically match species names in the `case.nml` to those in the mechanism. If you set `enable_reactions = .true.` (future feature), the solver can also discover and transport all species found in the mechanism automatically.

---

## 3. Best Practices
*   **Check for Typos**: The solver will now exit with a `fatal_error` if it encounters unknown variables in the namelist blocks.
*   **Peclet Number**: If using `"central"` convection, ensure your local Grid Peclet number ($|\mathbf{u}| \Delta x / \nu$) is $\le 2$ to avoid oscillations.
*   **Replicated Mesh**: Note that the mesh is currently replicated across all MPI ranks. For large meshes, ensure your nodes have sufficient memory.
