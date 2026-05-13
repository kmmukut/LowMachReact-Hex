title: Enthalpy/Species Coupling Convention

# Enthalpy/species coupling convention

## Decision

The solver uses **Option A**:

```text
h is the transported thermodynamic state.
T is recovered from h, Y, and p0.
```

When passive species transport updates the composition from `Y_old` to `Y_new`,
the solver must preserve the transported enthalpy field and recover temperature
from the new thermodynamic state:

```text
T_new = T(h_transported, Y_new, p0)
```

The solver must not preserve old temperature by recomputing enthalpy as:

```text
h_new = h(T_old, Y_new, p0)
```

That second form can numerically add or remove sensible enthalpy when composition
changes without an energy flux/source that accounts for it.

## Time-step convention

For the current constant-density, non-reacting enthalpy solver:

```text
1. Advance momentum/projection.
2. Advance passive species, if enabled.
3. Preserve the transported enthalpy h.
4. If Cantera thermo is enabled and transported species may have changed Y:
      sync T/cp/lambda/rho_thermo from h, current Y, and p0
      restore h so the thermo sync cannot alter the transported state
5. Advance h with advection, conduction, and qrad/rho.
6. If Cantera thermo is enabled:
      sync T/cp/lambda/rho_thermo from the updated h and current Y
      restore h so the thermo sync cannot alter the transported state
7. If Cantera thermo is disabled:
      recover T from h with the constant-cp relation.
```

For species-disabled Cantera thermo runs, step 4 may be skipped after initialization because composition has not changed and the previous post-energy sync already provides current `T`, `cp`, `lambda`, and `rho_thermo`.

```text
```

## Boundary convention

For fixed-temperature boundaries with species enabled, the boundary enthalpy is:

```text
h_b = h(T_b, Y_b, p0)
```

where `Y_b` comes from species boundary-condition evaluation. Dirichlet species
boundaries use `patch_Y`; non-Dirichlet species boundaries fall back to the
interior cell composition.

## Current limitations

This convention still omits the species-diffusion enthalpy correction:

```text
-div(sum_k h_k J_k)
```

Do not add reactions or variable-density low-Mach coupling until the passive
enthalpy/species tests pass.
