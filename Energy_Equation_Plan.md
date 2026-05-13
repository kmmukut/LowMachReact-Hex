# Enthalpy-Based Energy Equation Roadmap for Constant-Density Counterflow Validation

## Current status

The first passive sensible-enthalpy implementation is now in the codebase. This file remains as the development roadmap and validation plan, but several items that were originally future work are now implemented:

```text
- energy field storage for T, h, qrad, cp, lambda, and rho_thermo
- passive sensible-enthalpy transport
- Fourier conduction driven by grad(T)
- Cantera h(T,Y,p0) and T(h,Y,p0)
- combined Cantera thermo sync for T/cp/lambda/rho_thermo
- Option A enthalpy/species coupling: preserve h, recover T from h,Y,p0
- fixed-temperature boundary enthalpy using boundary composition
- qrad storage and diagnostics, with qrad initially zero
```

Still future work:

```text
- species-diffusion enthalpy correction: -div(sum_k h_k J_k)
- reactions and reaction heat release
- external radiation physics
- variable-density low-Mach projection
```

## Purpose

This note describes the staged plan for the enthalpy-based energy equation in the current finite-volume low-Mach codebase. The near-term target remains a **constant-density, non-reacting counterflow configuration** with accurate advection and diffusion of species and energy, validated against a non-reacting Cantera 1D opposed-flow/counterflow setup.

The longer-term target is a reacting, radiating low-Mach formulation where Cantera provides thermochemistry and transport, and an external radiation solver provides a volumetric source term `qrad`.

The recommended strategy is:

```text
Transport sensible enthalpy h as the energy variable.
Recover temperature T from h, Y, and thermodynamic pressure using Cantera.
Use T for transport properties, output, and validation.
Keep the first milestone constant-density and non-reacting.
Add radiation later as a volumetric enthalpy source term.
```

---

## Current Codebase Context

The current solver has the main passive scalar/energy pieces needed for this path:

1. **Finite-volume scalar transport pattern**
   - Species transport already provides a template for conservative face-based advection and diffusion.
   - The energy equation can follow this structure initially.

2. **Cantera bridge exists**
   - The bridge evaluates viscosity and mixture-averaged species diffusivities from Cantera.
   - The energy thermo path evaluates sensible enthalpy, recovers temperature from enthalpy, and refreshes `cp`, `lambda`, and diagnostic `rho_thermo`.

3. **Energy thermo state exists**
   - `mod_energy` owns `T`, `h`, `qrad`, `cp`, `lambda`, and `rho_thermo`.
   - The energy step uses a combined thermo sync where possible:
     ```text
     (T, cp, lambda, rho_thermo) = sync(h, Y, p0)
     ```

4. **Projection is currently constant-density/incompressible**
   - The flow solver uses a constant `rho` in the pressure equation and correction.
   - The current continuity target is effectively `div(u) = 0`.
   - Therefore, the first validation should not rely on thermal expansion or variable-density continuity.

The staged roadmap below respects this architecture while keeping the formulation compatible with future reacting and radiating flow.

---

## Why Use Enthalpy Instead of Temperature?

A temperature equation is simpler for the first implementation, but it becomes less natural once the solver includes:

- variable heat capacity,
- multi-species mixtures,
- species diffusion enthalpy flux,
- chemical heat release,
- radiation source terms,
- Cantera thermodynamic inversion,
- variable-density low-Mach coupling.

Sensible enthalpy is a better long-term transported variable because it is directly tied to energy conservation. Temperature should be treated as a thermodynamic state variable recovered from enthalpy and composition.

Use the relation:

```math
h = h(T, Y_1, Y_2, ..., Y_N, p_0)
```

and recover temperature from:

```math
T = T(h, Y_1, Y_2, ..., Y_N, p_0)
```

where `p0` is the thermodynamic pressure used for Cantera property evaluation. For the first constant-density target, `p0` can be fixed at the configured background pressure.

---

## Governing Variables

For the first enthalpy implementation, use cell-centered fields:

```text
u_c       velocity vector at cell center
p_c       hydrodynamic/projection pressure
Y_k,c     species mass fractions
h_c       sensible mixture enthalpy
T_c       temperature recovered from h, Y, p0
rho       constant density used by the flow solver
mu_c      dynamic viscosity
D_k,c     mixture-averaged species diffusivity
lambda_c  thermal conductivity
cp_c      mixture heat capacity at constant pressure
qrad_c    volumetric radiation source, initially zero
```

The first implementation should keep:

```text
rho = constant in the flow solver
p0 = constant thermodynamic pressure for Cantera
enable_reactions = false
qrad = 0
```

---

## Stage 1: Constant-Property Passive Enthalpy

Status: implemented.

### Goal

Validate the finite-volume advection and diffusion of the new energy scalar without introducing Cantera thermodynamic complexity.

### Model

Use constant density, constant heat capacity, and constant thermal conductivity:

```math
\rho \frac{D h}{D t} = \nabla \cdot (\lambda \nabla T)
```

with:

```math
h = c_p T
```

If `cp` is constant, this is equivalent to a passive temperature equation:

```math
\frac{D T}{D t} = \alpha \nabla^2 T
```

where:

```math
\alpha = \frac{\lambda}{\rho c_p}
```

### Finite-Volume Form

For each cell `c` with volume `V_c`:

```math
\rho V_c \frac{h_c^{n+1} - h_c^n}{\Delta t}
= - \sum_{f \in \partial c} \rho F_f h_f
  + \sum_{f \in \partial c} \lambda_f A_f \frac{T_{nb} - T_c}{d_f}
```

where:

```text
F_f   = u_f dot n_f A_f, volumetric face flux with owner-to-neighbor orientation
rho F_f = mass flux through the face
h_f   = upwind enthalpy at the face
a_f   = face area
d_f   = normal distance between cell centers, or cell center to boundary face
```

Because the current solver stores volumetric flux, use `rho * F_f` for enthalpy advection.

### Face Advection

Use upwind reconstruction for the first robust version:

```math
h_f =
\begin{cases}
h_c, & F_{out,c} \ge 0 \\
h_{nb}, & F_{out,c} < 0
\end{cases}
```

where `F_out,c` is the face flux oriented outward from cell `c`.

This should match the current species transport strategy and is the safest initial discretization.

### Face Diffusion

Use central diffusion in temperature, not enthalpy:

```math
q_{cond,f} = -\lambda_f \nabla T \cdot n_f
```

For a cell-centered finite-volume discretization:

```math
\left(\nabla T \cdot n\right)_f \approx \frac{T_{nb} - T_c}{d_f}
```

The contribution to cell `c` can be written as:

```math
\lambda_f A_f \frac{T_{nb} - T_c}{d_f}
```

Use arithmetic interpolation initially:

```math
\lambda_f = \frac{1}{2}(\lambda_c + \lambda_{nb})
```

For strong discontinuities in conductivity, harmonic interpolation may be better later.

### Boundary Conditions

Add temperature/enthalpy boundary conditions with at least:

```text
Dirichlet: fixed inlet/wall temperature
Neumann: zero normal temperature gradient
```

For a Dirichlet temperature boundary:

```math
T_b = T_{specified}
```

and the boundary enthalpy should be initialized consistently:

```math
h_b = c_p T_b
```

For a zero-gradient boundary:

```math
T_b = T_c, \quad h_b = h_c
```

### Validation Tests

Before counterflow, validate with:

1. **1D thermal diffusion between fixed-temperature boundaries**
   - Zero velocity.
   - Constant `lambda`, `rho`, and `cp`.
   - Steady solution should be linear in temperature.

2. **Passive scalar advection**
   - Uniform velocity.
   - Prescribed inlet temperature.
   - Check boundedness and numerical diffusion.

3. **Advection-diffusion balance**
   - Known Peclet-number behavior.
   - Verify grid refinement reduces error.

### Expected Outcome

At the end of Stage 1, the solver should transport `h`, recover `T = h / cp`, and produce stable temperature fields equivalent to a passive-temperature equation.

---

## Stage 2: Cantera Thermodynamics Without Reactions

Status: implemented for passive sensible enthalpy and diagnostic thermodynamic properties.

### Goal

Keep the flow constant-density and non-reacting, but use Cantera for thermodynamic and transport properties.

### Model

Transport sensible enthalpy:

```math
\rho \frac{D h}{D t} = \nabla \cdot (\lambda \nabla T)
```

but now:

```math
h = h(T, Y, p_0)
```

and:

```math
T = T(h, Y, p_0)
```

Cantera should provide:

```text
T from h, Y, p0
h from T, Y, p0
cp(T, Y, p0)
lambda(T, Y, p0)
mu(T, Y, p0)
D_k(T, Y, p0)
optional rho_cantera(T, Y, p0) for diagnostics only at this stage
```

### Important Restriction

During this stage, do **not** use Cantera density to modify the projection equation.

Use:

```text
rho_flow = constant solver density
rho_cantera = diagnostic or future-use thermodynamic density
```

This prevents the energy addition from forcing a premature rewrite of the pressure projection method.

### Thermodynamic Update Sequence

The current explicit sequence is:

```text
1. Refresh transport properties on transport_update_interval.
2. Advance momentum/projection.
3. Advance passive species, if enabled.
4. Preserve transported h.
5. If species may have changed Y and Cantera thermo is enabled:
      sync T/cp/lambda/rho_thermo from h,Y,p0
      restore h
6. Advance enthalpy using current lambda and T gradients.
7. If Cantera thermo is enabled:
      sync T/cp/lambda/rho_thermo from updated h,Y,p0
      restore h
8. Output h, T, and thermodynamic properties.
```

The sync operation is logically equivalent to recovering temperature and refreshing properties, but the implementation may combine these into one Cantera pass and may reuse cached dependent thermo state when `h`, `Y`, and `p0` are unchanged.

### Cantera State Inversion

The preferred conceptual operation is:

```text
set gas state using h, p0, Y
read back T
```

This makes enthalpy the primary transported variable and temperature the dependent thermodynamic variable.

For initialization and boundary conditions, the inverse operation is also needed:

```text
set gas state using T, p0, Y
read back h
```

### Numerical Detail: Diffusion Uses Temperature Gradient

Even when `h` is transported, the Fourier heat conduction term should use `grad(T)`:

```math
\nabla \cdot (\lambda \nabla T)
```

not simply:

```math
\nabla \cdot (\alpha \nabla h)
```

The latter is only equivalent when `cp` is constant and composition effects are negligible.

### Cache Requirement

Transport-property caches must depend on:

```text
T
p0
Y_k
```

Energy thermo-sync caches must depend on the transported thermodynamic state:

```text
h_sensible
p0
Y_k
```

A thermo-sync cache may reuse only dependent fields:

```text
T
cp
lambda
rho_thermo
```

It must not overwrite or reinterpret the transported enthalpy field.

### Validation Tests

1. **Constant-composition variable-cp heating**
   - Prescribe an initial enthalpy/temperature distribution.
   - Recover temperature using Cantera.
   - Check that `h(T)` and `T(h)` are mutually consistent.

2. **No-flow thermal diffusion with Cantera lambda and cp**
   - Compare with a high-resolution 1D reference.

3. **Species mixing without reactions**
   - Two streams with different composition.
   - Verify `sum(Y_k)` remains close to one and bounded.
   - Verify temperature changes are consistent with mixture enthalpy and composition.

### Expected Outcome

At the end of Stage 2, the solver should support Cantera-backed non-reacting enthalpy transport while still using a constant-density flow projection.

---

## Stage 3: Non-Reacting Counterflow Validation

### Goal

Validate the coupled advection-diffusion of species and enthalpy in an opposed-flow/counterflow configuration against Cantera's 1D counterflow setup with reactions disabled.

### Configuration

Use two opposed inlets:

```text
left/fuel-side inlet:
  velocity or mass flux equivalent
  T_fuel
  Y_fuel,k

right/oxidizer-side inlet:
  velocity or mass flux equivalent
  T_oxidizer
  Y_oxidizer,k
```

Use non-reacting Cantera chemistry:

```text
enable_species = true
enable_reactions = false
enable_energy = true
```

At this stage, keep:

```text
rho_flow = constant
p0 = constant
qrad = 0
```

### What to Compare

Compare centerline profiles after reaching steady state:

```text
T(x)
h(x)
Y_F(x)
Y_O2(x)
Y_N2(x)
mixture fraction Z(x)
scalar dissipation trend
stagnation-plane location
```

The most useful first comparisons are:

1. **Temperature profile**
   - Does the thermal mixing layer have the correct thickness and shape?

2. **Species profiles**
   - Are the opposed stream species transported with the same numerical behavior as temperature/enthalpy?

3. **Mixture fraction**
   - Does mixture fraction behave as a conserved scalar in the non-reacting case?

4. **Stagnation plane**
   - Does the stagnation point occur where expected based on inlet momenta/flow rates?

### Caveat About Exact Agreement

The current solver is constant-density and enforces approximately:

```math
\nabla \cdot u = 0
```

A full Cantera counterflow formulation may include variable density and low-Mach expansion effects. Therefore, the first comparison should be interpreted as validation of:

```text
finite-volume scalar transport
Cantera property coupling
enthalpy-temperature inversion
counterflow boundary setup
non-reacting mixing behavior
```

not yet as exact validation of full variable-density low-Mach counterflow physics.

### Numerical Accuracy Requirements

For counterflow validation, monitor:

```text
max and RMS divergence
sum(Y_k) error
minimum Y_k
minimum and maximum T
enthalpy boundedness
net boundary flux
centerline profile convergence
steady residuals of T, h, and Y_k
```

Recommended scalar convergence checks:

```math
R_h = \frac{\|h^{n+1} - h^n\|_2}{\max(\|h^n\|_2, \epsilon)}
```

```math
R_{Y_k} = \frac{\|Y_k^{n+1} - Y_k^n\|_2}{\max(\|Y_k^n\|_2, \epsilon)}
```

A steady counterflow validation should not rely only on time step count; it should use scalar residuals and profile convergence.

---

## Stage 4: Add Volumetric Radiation Source `qrad`

### Goal

Add radiation as a volumetric source term in the enthalpy equation.

The radiation solver will provide:

```text
qrad_c [W/m^3]
```

at each flow cell.

### Sign Convention

Define the sign convention explicitly:

```text
qrad > 0 means radiation adds energy to the gas.
qrad < 0 means radiation removes energy from the gas.
```

If the radiation model reports radiative loss as a positive number, insert it as:

```math
q_{energy} = -q_{rad,loss}
```

### Energy Equation With Radiation

For non-reacting flow:

```math
\rho \frac{D h}{D t}
= \nabla \cdot (\lambda \nabla T) + q_{rad}
```

Finite-volume form:

```math
\rho V_c \frac{h_c^{n+1} - h_c^n}{\Delta t}
= - \sum_f \rho F_f h_f
  + \sum_f \lambda_f A_f \frac{T_{nb} - T_c}{d_f}
  + q_{rad,c} V_c
```

### Numerical Treatment

The first implementation can treat radiation explicitly:

```math
h_c^{n+1} = h_c^n + \Delta t \left[ RHS_{adv+diff}^n + \frac{q_{rad,c}^n}{\rho} \right]
```

This is simple and consistent with an explicit scalar transport update.

If radiation becomes stiff, especially with strong gas cooling/heating, later options include:

```text
semi-implicit source treatment
operator splitting
source-term subcycling
coupled radiation-energy iteration
```

### Interface Recommendation

Design the energy module with a generic volumetric source array from the beginning:

```text
q_energy(:) [W/m^3]
```

Initially:

```text
q_energy = 0
```

Later:

```text
q_energy = qrad
```

Eventually:

```text
q_energy = qrad + qchem + other volumetric sources
```

This keeps the radiation coupling clean and avoids redesigning the energy solver interface.

---

## Stage 5: Add Reacting Enthalpy Physics

### Goal

After non-reacting scalar transport is validated, enable chemical source terms and heat release.

There are two common ways to formulate the reacting enthalpy equation.

### Option A: Transport Sensible Enthalpy With Chemical Source

Use:

```math
\rho \frac{D h_s}{D t}
= \nabla \cdot (\lambda \nabla T)
  - \nabla \cdot \left(\sum_k h_{s,k} j_k\right)
  + \dot{q}_{chem}
  + q_{rad}
```

where:

```text
h_s      sensible mixture enthalpy
h_s,k    sensible species enthalpy
j_k      diffusive species mass flux
qchem    chemical heat-release source
qrad     volumetric radiation source
```

### Option B: Transport Total Chemical Enthalpy

If total species enthalpies are used consistently, chemical heat release can be embedded through the species source terms and species enthalpies. This can be elegant, but it requires very careful consistency between:

```text
species equations
species production rates
species enthalpies
mixture enthalpy definition
Cantera reference states
```

### Recommendation

Do not start here.

First validate non-reacting enthalpy transport, then add reacting terms after the thermodynamic inversion and scalar numerics are reliable.

When reactions are added, prioritize consistency over speed:

```text
species source terms must match energy source terms
Cantera state must be updated from current T, Y, p0
enthalpy and temperature must remain thermodynamically consistent
```

---

## Stage 6: Variable-Density Low-Mach Coupling

### Goal

Move from passive constant-density energy/species transport to true low-Mach variable-density reacting flow.

This is a larger architectural change and should be delayed until the scalar transport, Cantera bridge, and counterflow validation are working.

### Why It Is Larger

The current projection enforces a divergence-free velocity field. Variable-density low-Mach flow generally requires a divergence constraint related to thermodynamic expansion and density changes.

Instead of:

```math
\nabla \cdot u = 0
```

one eventually needs a constraint of the form:

```math
\nabla \cdot u = S_{thermo}
```

where `S_thermo` depends on heat release, diffusion, composition change, pressure evolution, and possibly radiation.

### Required Changes

A true variable-density low-Mach upgrade would affect:

```text
pressure projection RHS
velocity correction coefficient
face flux definition
species transport fluxes
energy transport fluxes
CFL calculation
boundary mass balance
Cantera density coupling
output diagnostics
```

### Recommendation

Use Cantera density as diagnostic data first:

```text
rho_cantera(T,Y,p0)
```

but keep:

```text
rho_flow = constant
```

until the variable-density projection is deliberately designed and tested.

---

## Numerics Summary

### Transported Quantities

Use conservative finite-volume updates for:

```text
Y_k
h
```

Recover:

```text
T = T(h,Y,p0)
```

Evaluate:

```text
mu(T,Y,p0)
D_k(T,Y,p0)
lambda(T,Y,p0)
cp(T,Y,p0)
```

### Recommended Initial Time Integration

Use the same explicit style as the current scalar transport:

```text
Euler first
then optionally AB2 after validation
```

For a scalar `phi`:

```math
\phi^{n+1} = \phi^n + \Delta t \; RHS^n
```

For enthalpy:

```math
h^{n+1} = h^n + \Delta t
\left[
- \frac{1}{V_c} \sum_f F_f h_f
+ \frac{1}{\rho V_c} \sum_f \lambda_f A_f \frac{T_{nb} - T_c}{d_f}
+ \frac{q_{energy,c}}{\rho}
\right]
```

where `F_f` is the outward volumetric flux for cell `c`.

### Stability Limits

The explicit scalar update should respect both advective and diffusive limits.

Advective CFL:

```math
C_{adv,c} = \frac{\Delta t}{2 V_c} \sum_f |F_f|
```

Thermal diffusion limit estimate:

```math
C_{diff,T,c} = \Delta t \; \alpha_c \sum_f \frac{A_f}{V_c d_f}
```

where:

```math
\alpha_c = \frac{\lambda_c}{\rho c_{p,c}}
```

Species diffusion limit estimate:

```math
C_{diff,Y,c,k} = \Delta t \; D_{k,c} \sum_f \frac{A_f}{V_c d_f}
```

A conservative first target is:

```text
max(C_adv) < 0.5
max(C_diff_T) < 0.5
max(C_diff_Y) < 0.5
```

The exact limit depends on dimensionality, mesh quality, and discretization details, so these should be treated as conservative engineering limits, not mathematical guarantees.

### Boundedness

For early validation, enforce or monitor:

```text
T_min <= T <= T_max physically reasonable
0 <= Y_k <= 1
sum_k Y_k approx 1
h consistent with T and Y after Cantera inversion
```

If clipping is used, log it. Silent clipping can hide numerical problems.

### Species Normalization

Species transport should preserve:

```math
\sum_k Y_k = 1
```

In practice, numerical drift can occur. For non-reacting validation, monitor:

```math
E_Y = \max_c \left| \sum_k Y_{k,c} - 1 \right|
```

If renormalization is used, it should be done consistently and reported.

### Enthalpy-Temperature Consistency

After every enthalpy update:

```text
set Cantera state from h, Y, p0
read T
optionally read h_check from T, Y, p0
verify h_check approx h
```

Track:

```math
E_h = \max_c \frac{|h_{check,c} - h_c|}{\max(|h_c|, \epsilon)}
```

This is a critical diagnostic for the Cantera coupling.

---

## Recommended Output Fields

Add the following fields to VTU output for validation:

```text
T
h
rho_flow
rho_cantera
mu
lambda
cp
alpha
qrad
sum_Y
selected Y_k
```

For diagnostics CSV, add:

```text
min_T
max_T
mean_T
min_h
max_h
max_sumY_error
max_enthalpy_inversion_error
max_scalar_residual_T_or_h
max_species_residual
```

For counterflow validation, also output centerline/sample-line profiles in CSV:

```text
x, y, z, T, h, Y_F, Y_O2, Y_N2, Z, velocity components
```

This will make comparison against Cantera 1D much easier than relying only on visualization files.

---

## Recommended Development Order

Use this order to minimize risk:

1. Add energy field storage and initialization.
2. Add temperature boundary-condition parsing and evaluation.
3. Add constant-property enthalpy transport.
4. Add VTU and CSV output for `h` and `T`.
5. Validate pure thermal diffusion.
6. Couple Cantera `T <-> h` conversion. (Done)
7. Add Cantera `cp`, `lambda`, `mu`, and `D_k` updates from cellwise `T,Y,p0`. (Done for current passive modes)
8. Add combined Cantera thermo sync and conservative cache. (Done)
9. Validate no-flow Cantera thermal diffusion.
10. Run non-reacting opposed-flow/counterflow with species and enthalpy.
11. Add manufactured `qrad` source tests.
12. Couple external radiation source as volumetric `qrad`.
13. Only after this, design variable-density low-Mach projection.

---

## Definition of Done for the First Milestone

The first milestone is complete when the solver can run:

```text
constant-density
non-reacting
species enabled
enthalpy enabled
radiation disabled
Cantera thermodynamics enabled
opposed-flow/counterflow boundary conditions
```

and produce stable, converged profiles for:

```text
T(x)
h(x)
Y_k(x)
sum(Y_k)
```

with diagnostics showing:

```text
bounded species
bounded temperature
small sum(Y_k) error
small enthalpy-temperature inversion error
stable scalar residuals
reasonable agreement with non-reacting Cantera 1D reference trends
```

---

## Key Design Decisions

### Decision 1: Transport h, not T

Transporting enthalpy now avoids a later rewrite when adding radiation, reactions, and variable heat capacity.

### Decision 2: Diffuse T, not h

The conductive heat flux is proportional to `grad(T)`, so the diffusion term should use temperature gradients even though the transported variable is enthalpy.

### Decision 3: Keep density constant initially

This isolates scalar transport and Cantera thermodynamics from the harder problem of variable-density projection.

### Decision 4: Add qrad as a volumetric source

Because the radiation model will provide `qrad [W/m^3]`, the flow solver should consume it as a local source term in the enthalpy equation.

### Decision 5: Treat Cantera density as diagnostic first

Do not feed Cantera density into the pressure projection until the low-Mach divergence constraint is redesigned.

---

## Final Recommended Near-Term Equation Set

For the next implementation target, solve:

### Momentum and Projection

Current constant-density projection system:

```math
\nabla \cdot u = 0
```

using the existing pressure-correction method.

### Species

For each transported species:

```math
\frac{\partial Y_k}{\partial t} + \nabla \cdot (u Y_k)
= \nabla \cdot (D_k \nabla Y_k)
```

with no reactions initially.

### Enthalpy

```math
\rho \frac{\partial h}{\partial t}
+ \rho \nabla \cdot (u h)
= \nabla \cdot (\lambda \nabla T) + q_{energy}
```

where initially:

```math
q_{energy} = 0
```

and later:

```math
q_{energy} = q_{rad}
```

### Thermodynamic Closure

```math
T = T(h, Y, p_0)
```

using Cantera.

This gives a clean, practical path from passive temperature-like validation to full enthalpy-based reacting/radiating low-Mach simulation.
