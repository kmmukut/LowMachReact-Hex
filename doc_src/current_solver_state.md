title: Current Solver State

# Current Solver State: Energy, Species, and Cantera Thermo

## Summary

The current solver should be treated as a constant-density low-Mach / incompressible projection code with optional passive species transport and optional passive sensible-enthalpy transport.

The projection density remains the configured flow density:

```text
transport%rho = params%rho
```

Cantera thermodynamic density is stored separately as `energy%rho_thermo` and is diagnostic / future-use only. It must not be fed into the pressure projection until a separate variable-density low-Mach formulation is designed and validated.

## Main runtime ordering

The intended time-step ordering is:

```text
1. Refresh transport properties on transport_update_interval.
   - Cantera mu and D_k may use energy%T when energy is enabled.
   - Flow density remains params%rho.
2. Advance momentum / pressure projection.
3. Advance species transport, if enabled.
4. Advance sensible enthalpy transport, if enabled.
   - With Cantera thermo, sync T/cp/lambda/rho_thermo from h,Y,p0.
   - If transported species are enabled, perform a pre-flux thermo sync after species transport.
   - If species are disabled, the pre-flux thermo sync may be skipped because the previous post-flux sync is still valid.
   - Always perform the post-flux thermo sync after h is updated.
5. Write diagnostics and VTU/PVTU output on output steps.
```

## Current enabled physics

Supported in the current staged implementation:

```text
- constant flow/projection density
- passive species mass-fraction transport
- passive sensible-enthalpy transport
- Fourier conduction driven by grad(T)
- Cantera h(T,Y,p0), T(h,Y,p0), cp, lambda, and diagnostic rho_thermo
- qrad storage and diagnostics with qrad initially zero unless a future coupling fills it
```

Not yet implemented:

```text
- variable-density low-Mach divergence constraint
- reactions and reaction heat release
- species-diffusion enthalpy correction: -div(sum_k h_k J_k)
- external radiation physics
- radiation MPI/spectral decomposition
```

## Boundary-state rule

## Thermo update interval rule

`thermo_update_interval` is reserved for future optimization. The current supported value is:

```text
thermo_update_interval = 1
```

The solver may optimize internally using a combined thermo-sync call and conservative cache, but the logical thermodynamic state must remain synchronized every energy step.

For a fixed-temperature boundary with species enabled and Cantera thermo enabled, boundary enthalpy must use the boundary composition:

```text
h_b = h(T_b, Y_b, p0)
```

where `Y_b` comes from the species boundary condition, such as `patch_Y` on fixed-value inlet patches. This avoids using an interior composition or a default bath-gas composition for fuel/oxidizer inlet enthalpy.
