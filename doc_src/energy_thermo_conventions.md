title: Energy and Thermodynamic Conventions

# Energy and Thermodynamic Conventions

## Transported variable

The transported energy variable is mixture sensible enthalpy:

```text
h [J/kg]
```

Temperature is a dependent thermodynamic variable:

```text
T = T(h, Y, p0)
```

For the current constant-density implementation, `p0` is:

```text
p0 = params%background_press
```

Do not use projection pressure as thermodynamic pressure in the current stage.

## Sensible enthalpy reference

Cantera returns absolute mixture enthalpy including formation contributions. The current non-reacting energy model uses sensible enthalpy relative to `energy_reference_T` at the same composition:

```text
h_sens(T,Y,p0) = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
```

Temperature recovery performs the inverse by adding the same-composition reference enthalpy before HP inversion:

```text
h_abs_target = h_sens + h_abs(T_ref,Y,p0)
T = T(h_abs_target,p0,Y)
```

This prevents artificial heat release when passive species composition changes without reactions.

## Heat conduction

## Thermo sync and caching

When Cantera thermo is enabled, the preferred energy-step operation is a combined sync:

```text
(T, cp, lambda, rho_thermo) = thermo_sync(h, Y, p0)
```

This is equivalent to recovering temperature from enthalpy and then refreshing properties, but it is performed in one Cantera cell loop where possible.

The combined sync must preserve the transported enthalpy:

```text
h_after_sync = h_before_sync
```

A conservative cache may reuse the dependent thermo state only when the cache key is unchanged:

```text
h_sens, Y, p0
```

Do not use a cache keyed only on `T` or only on `Y` for energy thermodynamics.

Even though `h` is transported, thermal diffusion is driven by temperature gradient:

```text
conduction uses grad(T), not grad(h)
```

The current finite-volume update should be interpreted as:

```text
rho * D h / D t = div(lambda grad T) + qrad
```

where `rho = params%rho`.

## Radiation source convention

Use this sign convention everywhere:

```text
qrad > 0: radiation adds energy to the gas
qrad < 0: radiation removes energy from the gas
```

If a radiation model reports positive radiative loss, convert it in the coupling layer:

```text
qrad = -q_loss
```

## Units

```text
h       J/kg
T       K
rho     kg/m^3
cp      J/kg/K
lambda  W/m/K
mu      Pa s
nu      m^2/s
D_k     m^2/s
qrad    W/m^3
p0      Pa
```
