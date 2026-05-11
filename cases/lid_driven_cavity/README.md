# Lid-Driven Cavity

Generate the mesh from an activated `react_env`:

```bash
./make_mesh.sh
```

Run from the project root:

```bash
make
mpirun -np 1 ./hex_lowmach_fv cases/lid_driven_cavity/case.nml
```

The top `ymax` wall moves with `u = 1 m/s`. All other walls are stationary.

