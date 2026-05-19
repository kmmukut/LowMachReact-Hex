title: Energy, Thermodynamic, and Enthalpy/Species Coupling Conventions

# Energy, Thermodynamic, and Enthalpy/Species Coupling Conventions

This document defines the energy-variable, thermodynamic-state, and
enthalpy/species coupling conventions used by LowMachReact-Hex.

It merges the previous energy/thermodynamic convention note with the
enthalpy/species coupling decision.  The merged document reflects the current
implementation:

```text
stable baseline:
  constant-density sensible-enthalpy transport

active development layer:
  guarded non-reacting variable-density low-Mach sensible-enthalpy transport
  using Cantera thermo density and properties

future extensions:
  reaction source terms, reaction heat release, and radiation coupling through qrad
```

The solver transports sensible enthalpy.  Temperature, heat capacity, thermal
conductivity, and thermodynamic density are dependent properties recovered from
the transported state, composition, and thermodynamic pressure.

---

## 1. Primary decision: Option A

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

The solver must not preserve old temperature by rebuilding enthalpy as:

```text
h_new = h(T_old, Y_new, p0)
```

That second form can numerically add or remove sensible enthalpy when composition
changes without an energy flux or source term that accounts for it.

This convention applies to both:

```text
constant-density energy transport:
  h is the transported state

variable-density energy transport:
  rho*h is updated conservatively, then h remains the thermodynamic state used
  for Cantera temperature recovery
```

---

## 2. Transported energy variable

The transported thermodynamic variable is mixture sensible enthalpy:

```text
h [J/kg]
```

Temperature is a dependent thermodynamic variable:

```text
T = T(h, Y, p0)
```

where:

```text
Y   = transported species mass fractions or default/bath composition
p0  = thermodynamic/background pressure
```

The energy equation should not treat temperature as the primary transported
state when Cantera thermo is enabled.  Temperature is recovered from enthalpy and
composition.

---

## 3. Thermodynamic pressure

The thermodynamic pressure used for Cantera state calls is:

```text
p0 = params%background_press
```

This is distinct from the projection pressure field.

```text
background_press:
  thermodynamic pressure used in Cantera TPY/HP state evaluations

projection pressure:
  hydrodynamic pressure-like field used by the pressure projection
```

Do not use the projection pressure as the Cantera thermodynamic pressure unless
the low-Mach formulation is deliberately changed and revalidated.

High-pressure or real-gas cases should set `background_press` consistently with
the selected Cantera phase and mechanism.

---

## 4. Sensible enthalpy reference convention

Cantera absolute mixture enthalpy includes formation contributions.  The current
non-reacting energy model uses sensible enthalpy relative to
`energy_reference_T` at the same composition:

```text
h_sens(T,Y,p0) = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
```

where:

```text
T_ref = params%energy_reference_T
```

Temperature recovery performs the inverse operation by adding the
same-composition reference enthalpy before Cantera HP inversion:

```text
h_abs_target = h_sens + h_abs(T_ref,Y,p0)
T = T(h_abs_target,Y,p0)
```

This convention prevents artificial heat release when passive species
composition changes without reactions.

The same `T_ref` convention must be used consistently in:

```text
cantera_update_thermo_c
cantera_recover_temperature_from_h_c
cantera_recover_temperature_and_update_thermo_c
cantera_species_sensible_enthalpies_c
```

and in all matching Fortran `bind(c)` interfaces and call sites.

---

## 5. Composition convention

Thermodynamic properties are composition dependent.

For cell-centered thermo updates:

```text
if species are enabled:
  Y = species%Y(:, cell)

if species are disabled:
  Y = configured default/bath composition
```

For fixed-temperature and fixed-composition inlets:

```text
T_b = patch_T
Y_b = patch_Y
h_b = h(T_b,Y_b,p0)
rho_b = rho(T_b,Y_b,p0)
```

Dirichlet species boundaries use the configured boundary composition.  Non-
Dirichlet species boundaries use the interior composition for the boundary
state.

Do not evaluate fuel/oxidizer inlet enthalpy or density using an interior
composition or a default bath composition unless that is the intended boundary
state.

This boundary-state rule is important for:

```text
- fixed-temperature species inlets
- variable-density inlet density
- future mass-flow boundary conditions
- real-gas/high-pressure inlet states
- species-enthalpy diffusion consistency
```

---

## 6. Thermo sync

When Cantera thermo is enabled, the preferred energy-step operation is the
combined thermo sync:

```text
(T, cp, lambda, rho_thermo) = thermo_sync(h, Y, p0)
```

This is equivalent to:

```text
1. recover T from h,Y,p0
2. refresh cp, lambda, rho_thermo from T,Y,p0
```

but it is performed in one Cantera cell loop where possible.

The combined sync must preserve the transported enthalpy:

```text
h_after_sync = h_before_sync
```

The sync updates only dependent thermodynamic properties:

```text
T
cp
lambda
rho_thermo
```

It must not rebuild `h` from the recovered temperature after species changes.

Implementation rule:

```text
If a thermo-sync call could alter h through an older path or temporary state,
restore the transported h after the sync.  The transported h array remains the
source of truth.
```

---

## 7. Time-step convention

### 7.1 Shared ordering

The current high-level time-step ordering is:

```text
1. Advance momentum/projection.
2. Advance passive species, if enabled.
3. Advance sensible enthalpy, if enabled.
4. If variable-density mode is enabled:
     sync active density from thermo
     update the low-Mach divergence source for the next projection
```

The energy routine owns the detailed thermo-sync sequence inside step 3.

### 7.2 Constant-density Cantera-thermo energy ordering

For the constant-density, non-reacting enthalpy solver with Cantera thermo:

```text
1. Species transport may update Y.
2. Preserve the transported enthalpy h.
3. If transported species are enabled and Y may have changed:
      sync T/cp/lambda/rho_thermo from h, current Y, and p0
      preserve/restore h so the thermo sync cannot alter the transported state
4. Advance h with:
      advection
      conduction using grad(T)
      optional species-enthalpy diffusion
      qrad/rho
5. Sync T/cp/lambda/rho_thermo from updated h, current Y, and p0.
6. Preserve/restore h so the thermo sync cannot alter the transported state.
```

For species-disabled Cantera thermo runs, the pre-flux thermo sync may be skipped
after initialization because composition has not changed and the previous
post-flux sync already provides current `T`, `cp`, `lambda`, and `rho_thermo`.

### 7.3 Variable-density Cantera-thermo energy ordering

In variable-density low-Mach mode, the same Option A thermodynamic convention
applies, but the conservative update is in terms of `rho*h`.

The variable-density branch should be interpreted as:

```text
1. Use active density history and current fields to form the conservative update.
2. Advance the conservative energy state rho*h.
3. Recover h from the updated conservative state and active density.
4. Sync T/cp/lambda/rho_thermo from h, current Y, and p0.
5. Sync active flow density:
      transport%rho <- energy%rho_thermo
6. Update the low-Mach source for the next projection.
```

The energy thermo sync still must preserve the transported thermodynamic state
`h`.  The density update changes the active flow density for later flow,
species, and energy operations; it should not be interpreted as a license to
rebuild enthalpy from temperature.

---

## 8. Cache dependency rules

Cantera calls may be cached, but cache keys must include the actual independent
thermodynamic state.

For transport properties:

```text
cache key: T, p0, Y_1...Y_N
outputs:   mu, D_k
```

For combined energy thermo sync:

```text
cache key: h_sens, p0, Y_1...Y_N
outputs:   T, cp, lambda, rho_thermo
```

For species sensible enthalpies:

```text
cache key: T, p0, Y_1...Y_N
outputs:   h_k(T,Y,p0) - h_k(T_ref,Y,p0)
```

Do not cache energy thermodynamics only on `T` or only on `Y`.  The transported
state is `h_sens`, and the composition and thermodynamic pressure are part of
the state.

A cache hit may reuse dependent properties, but it must not reinterpret or
overwrite the transported enthalpy field.

---

## 9. Heat conduction convention

Even though `h` is transported, thermal diffusion is driven by the temperature
gradient:

```text
conduction uses grad(T), not grad(h)
```

Constant-density interpretation:

```text
rho * D h / D t = div(lambda grad T) + qrad + species_enthalpy_terms
```

where:

```text
rho = params%rho
```

Variable-density conservative interpretation:

```text
\partial(\rho h)/\partial t + \nabla\cdot(\rho u h)
  = \nabla\cdot(\lambda \nabla T) + qrad + species_enthalpy_terms
```

The finite-volume conduction contribution is a heat flux driven by temperature
difference across faces:

```text
q_cond ~ -lambda grad(T)
```

For diagnostics, the outward conductive boundary flux is positive when energy
leaves the domain.

---

## 10. Species-enthalpy diffusion convention

When enabled, the species-enthalpy diffusion correction accounts for enthalpy
transported by diffusive species fluxes:

```text
-div(sum_k h_k J_k)
```

where `h_k` is species sensible enthalpy relative to the same `T_ref`:

```text
h_k_sens = h_k(T,Y,p0) - h_k(T_ref,Y,p0)
```

This term should be evaluated consistently with the species diffusion fluxes
used by the species transport operator.

In variable-density mode, the species diffusion and species-enthalpy diffusion
paths should be density-weighted consistently with conservative `rho*Y` and
`rho*h` transport.

Recommended validation toggles:

```fortran
enable_species_enthalpy_diffusion = .false.
enable_species_enthalpy_diffusion = .true.
```

Use both to distinguish base enthalpy-transport behavior from the additional
species-enthalpy correction.

---

## 11. Radiation source convention

The volumetric radiation/source field is:

```text
qrad [W/m^3]
```

Use this sign convention everywhere:

```text
qrad > 0: adds energy to the gas
qrad < 0: removes energy from the gas
```

If a radiation model reports positive radiative loss, convert it in the coupling
layer:

```text
qrad = -q_loss
```

The current solver provides the `qrad` storage and energy-equation hook.  A full
radiation solver is future work.  Manufactured `qrad` source tests should be
used before coupling a physical radiation model.

---

## 12. Density semantics

Cantera thermo sync returns:

```text
rho_thermo = rho(T,Y,p0)
```

How `rho_thermo` is used depends on solver mode.

### 12.1 Constant-density mode

```text
enable_variable_density = .false.
```

Active flow density remains:

```text
transport%rho = params%rho
```

Cantera density is diagnostic:

```text
energy%rho_thermo
```

The pressure projection, momentum update, mass flux, and global mass diagnostic
use `transport%rho`, not `rho_thermo`.

### 12.2 Variable-density low-Mach mode

```text
enable_variable_density = .true.
density_eos = "cantera"
```

Active flow density is synchronized from thermo:

```text
transport%rho <- energy%rho_thermo
```

The variable-density low-Mach path then uses this density for:

```text
- mass flux
- conservative species transport
- conservative enthalpy transport
- variable-coefficient projection
- low-Mach source construction
```

This path is active but should remain guarded and validation-driven until the
broader EOS/pressure/boundary-condition validation matrix is complete.

---

## 13. Variable-density low-Mach source convention

In variable-density low-Mach mode, the projection target is not raw
incompressibility.  The target is:

```text
div(u) = S
```

The current conservative source form is:

```text
S = (rho_old - rho)/(rho*dt) - (u.grad(rho))/rho
```

The first term captures local density change.  The second term accounts for
advected density gradients so that conservative continuity is targeted:

```text
\partial\rho/\partial t + \nabla\cdot(\rho u) = 0
```

Projection validation should use the source time level actually consumed by the
pressure solve:

```text
divu_minus_S_projection_*
```

Current-source residuals are useful for diagnosing source evolution after the
energy/thermo update, but they are not the primary projection pass/fail metric.

---

## 14. Energy budget diagnostics

The energy diagnostics distinguish several concepts that should not be mixed.

### 14.1 Direct energy-update closure

The direct update closure checks whether the finite-volume update itself closes:

```text
rho_new h_new = rho_old h_old + dt * RHS
```

Primary metric:

```text
relative_last_energy_update_balance_defect
```

This should be near roundoff for a clean operator update.

### 14.2 Operator-consistent rho-h budget

The operator-consistent budget uses the density/time level actually used by the
energy operator.

Primary metric:

```text
rel_operator_recon_defect
```

### 14.3 Output-state reconciliation

The output-state snapshot may use a density/time level that differs from the
operator density/time level.  Reconciliation diagnostics explain this
bookkeeping difference.

Primary metric:

```text
rel_output_recon_defect
```

Older unreconciled endpoint budget columns are useful context, but should not be
used alone as conservation failure metrics.

---

## 15. Units

```text
h                 J/kg
h_k               J/kg_k
T                 K
rho               kg/m^3
cp                J/kg/K
lambda            W/m/K
mu                Pa s
nu                m^2/s
D_k               m^2/s
qrad              W/m^3
p0                Pa
u                 m/s
face_flux         m^3/s
mass_flux         kg/s
rho*h             J/m^3
rho*h integrated  J
```

---

## 16. Developer rules

When modifying the energy or thermo path:

```text
- preserve the transported variable convention: h is transported, T is recovered
- keep T_ref usage consistent across all Cantera calls
- distinguish background_press from projection pressure
- use boundary Y and T for boundary thermodynamic states
- do not let thermo sync overwrite transported h
- update cache keys if any independent state variable changes
- preserve constant-density behavior unless deliberately changing it
- add/update diagnostics for variable-density changes
- validate species-enthalpy diffusion both on and off
- validate species-enabled and species-disabled Cantera thermo paths separately
```

---

## 17. Future development goals

Near-term:

```text
- centralize boundary thermodynamic state evaluation
- expose boundary rho_b and h_b for variable-density inlet handling
- add patch-wise energy and mass-flow diagnostics
- keep energy diagnostics aligned with operator time levels
```

Short-to-medium term:

```text
- add mass-flow inlet BC using rho_b = EOS(T_b,Y_b,p0)
- broaden EOS/pressure validation for energy transport
- improve scalar/enthalpy convection accuracy with bounded higher-order options
- add restart support for long transient energy/variable-density runs
```

Long-term:

```text
- add chemical source terms and heat release after non-reacting validation
- couple physical radiation through qrad after manufactured qrad validation
- add stronger support for high-pressure real-gas boundary states
- consider HDF5/XDMF-style scientific output for restart and energy budgets
```

---

## 18. Status summary

```text
Option A h/Y/T coupling:
  active and central to the solver

Sensible enthalpy transport:
  active

Cantera T(h,Y,p0), cp, lambda, rho_thermo:
  active

Species-enthalpy diffusion:
  optional active correction

qrad source hook:
  active storage and energy-equation hook; physical radiation solver future

Constant-density energy path:
  stable baseline

Variable-density rho*h path:
  active and diagnostics-rich, still under broader validation

Reaction heat release:
  future work
```
