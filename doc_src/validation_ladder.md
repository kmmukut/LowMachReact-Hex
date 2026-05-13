title: Validation Ladder

# Validation Ladder for Energy, Species, and Cantera Thermo

Run these in order after each energy/thermo patch.

## 1. Build checks

```text
make clean
make release
make debug
```

Expected result: no Fortran interface mismatches and no missing `.mod` dependencies.

## 2. Energy-disabled regression

Run baseline cases with:

```text
enable_energy = .false.
```

Expected result: flow/species diagnostics match the pre-energy baseline within normal floating-point tolerance.

## 3. Energy initialization-only smoke test

Use:

```text
enable_energy = .true.
nsteps = 0 or a minimal run
qrad = 0
```

Expected result:

```text
T initialized to initial_T
h initialized consistently from T
qrad is zero
energy diagnostics and VTU/PVTU arrays are written
```

## 4. Cantera h-to-T roundtrip

With Cantera thermo enabled, verify:

```text
T_input -> h_sens(T,Y,p0) -> T_recovered
```

Expected result: recovered `T` matches the input within a small tolerance.

## 5. Fixed-temperature boundary composition test

Use at least two fixed-temperature inlet patches with different species compositions.

Expected result: boundary enthalpy is evaluated from:

```text
h_b = h(T_b, Y_b, p0)
```

not from the interior composition and not from the default bath-gas composition.

## 6. Pure diffusion test

Use zero velocity and fixed hot/cold temperature boundaries.

Expected result: temperature evolves smoothly toward the expected conductive profile; conduction is driven by `grad(T)`.

## 7. Pure advection test

Use weak or zero conduction and advect a smooth enthalpy/temperature profile.

Expected result: scalar transport follows the flow direction with bounded numerical diffusion.

## 8. Species + enthalpy non-reacting test

Enable species and energy, disable reactions, and keep `qrad = 0`.

Expected result:

```text
Y_k remains bounded
sum(Y_k) remains controlled
T remains physical
h/T/Y fields vary smoothly
```

## 9. Non-reacting counterflow comparison

Compare centerline trends against a non-reacting Cantera 1D counterflow reference.

Expected result: qualitative trends should agree, but exact agreement is not expected until variable-density low-Mach coupling is implemented.

## 10. Future qrad manufactured tests

Before external radiation physics, add prescribed source modes:

```text
uniform heating
uniform cooling
Gaussian heating/cooling
localized source
```

Expected result: the domain energy budget changes according to the volume integral of `qrad`.

## 11. Profiling regression for thermo-sync optimization

Use a Cantera thermo case with profiling enabled:

```text
enable_profiling = .true.
nested_profiling = .true.
```

Expected result after the combined thermo-sync optimization:

```text
Energy_Cantera_PreSync
Energy_Cantera_PostSync
```

replace the older four-region pattern:

```text
Energy_Cantera_PreRecoverT
Energy_Cantera_PreRefresh
Energy_Cantera_PostRecoverT
Energy_Cantera_PostRefresh
```

For species-disabled Cantera thermo runs, the pre-sync may be absent because composition does not change before the energy step.

Confirm that diagnostics remain physically equivalent before comparing runtime.

