# Test Plan

The current tree starts with validation-oriented smoke tests:

1. Activate `react_env`.
2. Build:

   ```bash
   make clean
   make
   ```

3. Generate meshes:

   ```bash
   ./cases/lid_driven_cavity/make_mesh.sh
   ./cases/channel_flow/make_mesh.sh
   ```

4. Run serial smoke tests:

   ```bash
   mpirun -np 1 ./lowmach_react_hex cases/lid_driven_cavity/case.nml
   mpirun -np 1 ./lowmach_react_hex cases/channel_flow/case.nml
   ```

5. Run MPI consistency smoke tests:

   ```bash
   mpirun -np 2 ./lowmach_react_hex cases/lid_driven_cavity/case.nml
   mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
   ```

Useful first checks:

- `diagnostics.csv` exists.
- `.vtu` visualization files exist and open in ParaView.
- `max_divergence` decreases after projection.
- cavity `net_boundary_flux` stays near zero.
- channel output remains finite with periodic `x` and `z` patches.
