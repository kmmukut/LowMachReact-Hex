title: Cantera Interface Notes


# Cantera Interface Notes

This document describes the current solver-facing Cantera interface rules.  It
complements the C++ bridge source note by focusing on how the Fortran solver
should call and interpret Cantera thermo/transport data.

The current interface supports both the stable constant-density solver path and
the guarded experimental non-reacting variable-density low-Mach path.

---

## 1. Interface role

`src/cantera_interface.cpp` bridges the Fortran solver to Cantera for
thermodynamic and transport properties.

The interface is responsible for:

```text
- loading a Cantera mechanism and optional named phase
- mapping solver species names to Cantera species indices
- building full Cantera mass-fraction vectors from transported species
- filling missing composition with a bath/default species when needed
- evaluating viscosity and scalar species diffusivities
- evaluating sensible mixture enthalpy h_sens(T,Y,p0)
- recovering temperature from transported sensible enthalpy h_sens
- refreshing cp, lambda, and rho_thermo
- evaluating species sensible enthalpies for species-enthalpy diffusion
- reporting bridge cache statistics
```

The interface is **not** the owner of solver state.  It evaluates properties
from state supplied by Fortran.  Flow, species, energy, and transport arrays
remain Fortran-owned.

---

## 2. Solver state versus dependent properties

The transported energy state is sensible enthalpy:

```text
h_sens
```

Temperature is a dependent thermodynamic property recovered from:

```text
T = T(h_sens, Y, p0)
```

Cantera thermo sync updates dependent fields:

```text
T
cp
lambda
rho_thermo
```

It must not reinterpret or overwrite the transported enthalpy field.

In constant-density mode:

```text
transport%rho = params%rho
energy%rho_thermo = diagnostic thermodynamic density
```

In variable-density mode:

```text
transport%rho <- energy%rho_thermo
```

The decision to use `rho_thermo` as active density belongs to the Fortran
solver mode, not to the C++ bridge itself.

---

## 3. Cantera mechanism and phase selection

The interface supports both a mechanism file and an optional phase name:

```fortran
&fluid_input
  cantera_mech_file  = "path/to/mechanism.yaml"
  cantera_phase_name = "phase_name"
/
```

Blank phase name preserves Cantera's default/first-phase behavior.  A nonblank
phase name selects a named phase from the YAML file.

The solver-level selector:

```fortran
density_eos = "cantera"
```

means:

```text
use density returned by the selected Cantera phase
```

It does **not** mean one specific EOS.  The actual thermodynamic model is set by
the selected Cantera phase, for example:

```yaml
thermo: ideal-gas
```

```yaml
thermo: Peng-Robinson
```

```yaml
thermo: Redlich-Kwong
```

Similarly, the selected phase determines the Cantera transport model, for
example:

```yaml
transport: mixture-averaged
transport: multicomponent
transport: high-pressure
transport: high-pressure-Chung
```

For the current solver path, gas-mixture phases are the validated target.  Other
Cantera phase classes may require additional solver assumptions and validation.

---

## 4. Thermodynamic pressure convention

Cantera state calls use the thermodynamic/background pressure supplied by the
Fortran side, normally:

```text
p0 = params%background_press
```

This pressure is separate from the projection pressure field.

```text
background_press:
  thermodynamic pressure used in Cantera TPY/HP state evaluations

projection pressure:
  hydrodynamic pressure-like variable used by the pressure projection
```

High-pressure EOS cases should set `background_press` consistently with the
chosen Cantera phase and mechanism.  Fixed pressure boundaries in the projection
solve do not replace the thermodynamic pressure used in Cantera property calls.

---

## 5. Sensible enthalpy convention

The interface uses sensible enthalpy relative to a reference temperature:

```text
h_sens(T,Y,p0) = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
```

The temperature recovery path reconstructs the absolute target enthalpy before
using Cantera's HP state solve:

```text
h_abs_target = h_sens + h_abs(T_ref,Y,p0)
setState_HP(h_abs_target, p0)
```

This convention avoids artificial heat release in non-reacting species mixing.
Composition changes alter the relationship between `h_sens` and `T`, but the
transported `h_sens` field remains the conserved/transported state.

The same `T_ref` convention must be used consistently by:

```text
cantera_update_thermo_c
cantera_recover_temperature_from_h_c
cantera_recover_temperature_and_update_thermo_c
cantera_species_sensible_enthalpies_c
```

and by the matching Fortran `bind(c)` interfaces and call sites.

---

## 6. Combined thermo-sync path

The energy step should prefer the combined thermo-sync bridge whenever it needs
both temperature recovery and refreshed thermodynamic properties:

```text
cantera_recover_temperature_and_update_thermo_c(
    h_sens,
    p0,
    Y,
    T_out,
    cp_out,
    lambda_out,
    rho_thermo_out,
    T_ref,
    species_names_flat,
    name_len
)
```

This call performs one cell loop and returns:

```text
T
cp
lambda
rho_thermo
```

It replaces the older two-pass pattern:

```text
recover_temperature_from_h
update_thermo_from_temperature
```

The older separate calls remain useful for initialization and point evaluations,
but the normal energy-step path should avoid doing a full HP inversion followed
by a full TPY property update when one synchronized pass is sufficient.

---

## 7. Composition rules

### 7.1 Cell thermo updates

For cell-centered thermo updates:

```text
if species are enabled:
  Y = species%Y(:, cell)

if species are disabled:
  Y = default/bath composition
```

The bridge maps solver species names to Cantera indices by name and builds a
full Cantera mass-fraction vector.

### 7.2 Missing species and bath composition

The bridge can fill missing composition with a bath species, currently preferring
`N2` when it exists.  This supports no-species thermo paths and reduced
transported species sets.

The solver should still keep transported species names and Cantera mechanism
species names consistent.  Missing names should be treated as a configuration
problem unless a deliberate reduced-species strategy is being used.

### 7.3 Boundary thermo states

For fixed-temperature/fixed-composition inlets, property evaluation should use
the boundary composition:

```text
Y_b = boundary species composition from patch_Y
T_b = boundary temperature from patch_T
```

Then:

```text
h_b   = h(T_b, Y_b, p0)
rho_b = rho(T_b, Y_b, p0)
```

This matters because both enthalpy and density are composition dependent.
Using an interior composition or a default bath composition for a fuel/oxidizer
inlet can give the wrong boundary state.

The current bridge can evaluate these properties if the Fortran side supplies
the correct boundary state.  A cleaner boundary-state helper API is a future
interface improvement.

---

## 8. Density rule

The bridge returns:

```text
rho_thermo = Cantera density for the supplied T/h, Y, and p0
```

How that density is used depends on solver mode.

### 8.1 Constant-density mode

```text
enable_variable_density = .false.
```

Active flow density remains:

```text
transport%rho = params%rho
```

Cantera density is:

```text
energy%rho_thermo
```

and should be treated as diagnostic/output/validation data.

### 8.2 Variable-density mode

```text
enable_variable_density = .true.
density_eos = "cantera"
```

Active density is synchronized from thermo:

```text
transport%rho <- energy%rho_thermo
```

The low-Mach projection and conservative species/energy branches then use this
active density.  This mode should be kept under guarded non-reacting validation
rules until the broader validation matrix is complete.

---

## 9. Transport-property rule

Transport-property calls must depend on the actual thermodynamic state:

```text
T
p0
Y_1 ... Y_N
```

Transport properties must not be cached only on composition once temperature or
pressure can evolve.

Current bridge outputs include:

```text
mu
D_k
lambda
```

Important interpretation:

```text
The FV species operator currently consumes scalar D_k values.  This naturally
matches mixture-averaged transport.  A Cantera phase may use a multicomponent
transport model, but the current FV operator does not yet assemble a full
coupled Stefan-Maxwell multicomponent diffusion flux.
```

For validation of the current species operator, mixture-averaged transport is
the closest conceptual match.  Multicomponent and high-pressure transport phases
remain useful for future work and selected property/EOS tests if the queried
properties are supported.

---

## 10. Cache dependency rules

Cantera calls are expensive, so the bridge includes conservative local caches.
The cache rules are part of the interface contract.

### 10.1 Transport cache

Transport-property cache key must depend on:

```text
T
p0
Y_1 ... Y_N
```

Cached outputs may include:

```text
mu
D_k
```

### 10.2 Combined thermo-sync cache

Combined thermo-sync cache key must depend on:

```text
h_sens
p0
Y_1 ... Y_N
```

Cached outputs may include:

```text
T
cp
lambda
rho_thermo
```

The cache must not overwrite or reinterpret the transported `h_sens` field.

### 10.3 Species-enthalpy diffusion cache

Species sensible enthalpy cache key must depend on:

```text
T
p0
Y_1 ... Y_N
```

Cached outputs may include:

```text
h_k(T,Y,p0) - h_k(T_ref,Y,p0)
```

This cache is only for species-enthalpy diffusion.  It must not overwrite or
reinterpret the transported mixture enthalpy field.

---

## 11. Species-enthalpy diffusion interface

The optional energy correction uses species sensible enthalpies:

```text
h_k_sens = h_k(T,Y,p0) - h_k(T_ref,Y,p0)
```

These are used in the finite-volume correction:

```text
-div(sum_k h_k J_k)
```

In variable-density mode, the species diffusion flux and enthalpy correction
should be density weighted consistently with the conservative `rho*Y` and
`rho*h` updates.

Validation should include runs with:

```text
enable_species_enthalpy_diffusion = .false.
enable_species_enthalpy_diffusion = .true.
```

to separate base enthalpy transport errors from species-enthalpy diffusion
effects.

---

## 12. Variable-density low-Mach interface contract

For the current non-reacting variable-density mode, the expected interface
contract is:

```text
1. Fortran transports species and enthalpy.
2. Cantera recovers T and updates cp/lambda/rho_thermo.
3. Fortran synchronizes active density:
      transport%rho <- energy%rho_thermo
4. Fortran updates low-Mach source:
      S = (rho_old - rho)/(rho*dt) - (u.grad(rho))/rho
5. The next projection enforces:
      div(u) = S_projection
```

The Cantera interface supplies thermodynamic density; the projection and
continuity consistency remain Fortran-side responsibilities.

---

## 13. Boundary-state interface gap

A recurring future need is a clean boundary-state property evaluation path.

Desired boundary helper semantics:

```text
given:
  T_b
  Y_b
  p0

return:
  h_b
  rho_b
  cp_b
  lambda_b
  mu_b
  D_k,b if needed
```

Use cases:

```text
- fixed-temperature species inlets
- density-aware variable-density inlet fluxes
- mass-flow boundary conditions
- real-gas/high-pressure inlet states
- enthalpy boundary consistency
```

Currently, these properties can be evaluated using existing bulk/point Cantera
paths if the Fortran side constructs the correct boundary arrays.  A dedicated
helper would make the interface clearer and reduce duplicated boundary logic.

---

## 14. Error and validation behavior

Cantera failures should be treated as configuration or state failures unless a
specific recovery strategy is implemented.

Common failure causes:

```text
- missing mechanism file
- wrong phase name
- unsupported transport model for a selected phase
- solver species not present in the mechanism
- invalid mass fractions
- HP recovery failure from out-of-range enthalpy/composition
- high-pressure EOS state outside the model's valid range
```

The current bridge exits on Cantera exceptions.  A future improvement is to
route these through the solver's MPI-aware fatal-error path so all ranks report
and terminate consistently.

---

## 15. Current limitations

The current interface does not yet provide:

```text
- reaction-rate source terms
- heat-release source terms
- full multicomponent diffusion flux matrices
- Soret or Dufour flux coupling
- dedicated boundary-state property helper
- multiple simultaneous active Cantera phases in one process
- thread-safe independent Cantera contexts
- nonfatal fallback for failed HP recovery
```

These limitations are acceptable for the current non-reacting validation track
if they remain explicit.

---

## 16. Future development assessment

### 16.1 Near term

```text
- remove duplicate/stale wording in older interface notes
- document named-phase and EOS behavior consistently
- add startup reporting for selected Cantera thermo and transport model
- add simple mechanism/phase probe tests
- route bridge errors through MPI-aware solver fatal handling
```

### 16.2 Short-to-medium term

```text
- add a dedicated boundary-state property evaluation helper
- use boundary rho(T_b,Y_b,p0) for variable-density inlet/outlet logic
- support mass-flow inlet BCs through boundary density
- expand validation across ideal-gas, Peng-Robinson, Redlich-Kwong, and
  pressure/transport variants
```

### 16.3 Medium-to-long term

```text
- add reaction source terms and heat release after non-reacting validation is robust
- add full multicomponent FV diffusion support only if the Fortran species
  operator is upgraded accordingly
- consider bridge context objects if multiple simultaneous mechanisms/phases
  become necessary
- consider solver-side ideal-gas EOS only if it provides clear value beyond Cantera
```

---

## 17. Status summary

```text
Cantera thermo and transport bridge:
  active

Named phase selection:
  active

Sensible enthalpy convention:
  active and central to energy transport

Cantera density:
  diagnostic in constant-density mode
  active in guarded variable-density mode

Species-enthalpy diffusion support:
  active as optional correction

Reaction and heat-release support:
  future work

Boundary density/enthalpy helper:
  recommended future interface improvement

Full multicomponent FV diffusion:
  future architectural extension, not current behavior
```