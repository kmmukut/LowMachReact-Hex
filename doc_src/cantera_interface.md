title: Cantera Interface Notes

# Cantera Interface Notes

## C++ bridge responsibilities

`src/cantera_interface.cpp` bridges the Fortran solver to Cantera for transport and thermodynamic properties.

The current staged responsibilities are:

```text
- map solver species names to Cantera species indices
- build a full Cantera mass-fraction vector for each cell
- fill missing composition with the configured bath/default species when needed
- compute mu and D_k for transport updates
- compute h_sens, cp, lambda, and rho_thermo for energy/thermo updates
- recover T from h_sens, Y, and p0
- recover T and refresh cp/lambda/rho_thermo in one combined thermo-sync pass for the energy step
- recover T and refresh cp/lambda/rho_thermo in one combined thermo-sync pass for the energy step
```

## Sensible enthalpy convention

The thermo update should compute:

```text
h_sens = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
```

The recovery path should compute:

```text
h_abs_target = h_sens + h_abs(T_ref,Y,p0)
setState_HP(h_abs_target, p0)
```

The `T_ref` argument must appear consistently in:

```text
cantera_update_thermo_c(..., rho_thermo_out, T_ref, species_names_flat, name_len)
cantera_recover_temperature_from_h_c(..., T_out, T_ref, species_names_flat, name_len)
```

and in the matching Fortran `bind(c)` interfaces and call sites.


## Combined thermo-sync path

The energy step should prefer the combined thermo-sync bridge when it needs both temperature recovery and refreshed properties:

```text
cantera_recover_temperature_and_update_thermo_c(
    h_sens, p0, Y,
    T_out, cp_out, lambda_out, rho_thermo_out,
    T_ref, species_names_flat, name_len
)
```

This call preserves the transported `h_sens` array and updates only dependent thermodynamic state:

```text
T
cp
lambda
rho_thermo
```

It replaces the older two-pass energy-step pattern:

```text
recover_temperature_from_h
update_thermo_from_temperature
```

The older separate calls remain useful for initialization and point evaluations, but the energy transport path should avoid doing a full HP inversion pass followed by a second full TPY property pass when one synchronized pass is sufficient.

## Composition rule

For cell thermo updates:

```text
Y = species%Y(:,cell) when species is enabled and present
Y = default inert/bath composition when species is disabled
```

For fixed-temperature species inlets:

```text
Y = boundary species composition from the boundary condition
```

This matters because `h(T,Y,p0)` is composition dependent.

## Density rule

`rho_thermo` is diagnostic only in the current solver. The projection, momentum update, and global mass diagnostic must continue using the constant flow density until variable-density low-Mach coupling is implemented deliberately.

## Cache dependency rule

Any Cantera cache for transport properties must depend on:

```text
T, p0, Y_1...Y_N
```

Do not cache transport properties only on composition once temperature evolves.

The combined energy thermo-sync cache must depend on the transported thermodynamic state:

```text
h_sens, p0, Y_1...Y_N
```

A cached thermo-sync result may reuse only the dependent fields:

```text
T, cp, lambda, rho_thermo
```

It must not overwrite or reinterpret the transported enthalpy field.
