# Rectangular Channel Flow

Generate the mesh from an activated `react_env`:

```bash
./make_mesh.sh
```

Run from the project root:

```bash
make
mpirun -np 1 ./hex_lowmach_fv cases/channel_flow/case.nml
```

The case uses periodic streamwise and spanwise patches, no-slip `y` walls,
and a constant streamwise body force.

