title: Cantera C++ Source Bridge


# Cantera C++ Source Bridge

The Cantera C++ bridge is implemented in:

```text
src/cantera_interface.cpp
```

This file is documented separately because FORD is primarily used for Fortran
source documentation, while the Cantera bridge is C++ code exposed to Fortran
through a C ABI.

The bridge is a process-local interface between the Fortran solver and Cantera.
Each MPI rank owns its own Cantera `Solution`, `ThermoPhase`, `Transport`
manager, and bridge-local caches.

---

## 1. Role in the solver

`cantera_interface.cpp` provides the C-callable wrappers used by the Fortran
solver for thermodynamics and transport.

The bridge currently supports:

```text
- mechanism and phase initialization
- species-count and species-name queries
- solver-species to Cantera-species mapping
- construction of full Cantera mass-fraction vectors
- dynamic viscosity evaluation
- scalar species diffusivity evaluation
- sensible mixture enthalpy evaluation
- temperature recovery from transported sensible enthalpy
- combined temperature recovery and thermo-property synchronization
- species sensible enthalpy evaluation for species-enthalpy diffusion
- cache-effectiveness diagnostics
```

The bridge does **not** own the transported solution state.  Solver state remains
in Fortran arrays.  The C++ bridge is called to evaluate dependent properties.

Core thermodynamic sync:

```text
given:  h, Y, p0
return: T, cp, lambda, rho_thermo
```

The transported energy variable remains sensible enthalpy `h`.  The bridge
updates dependent thermodynamic properties only.

---

## 2. Process-local Cantera objects

The bridge stores process-local global pointers:

```cpp
static std::shared_ptr<Cantera::Solution> sol;
static std::shared_ptr<Cantera::ThermoPhase> gas;
static std::shared_ptr<Cantera::Transport> trans;
```

Interpretation:

```text
- every MPI rank initializes its own Cantera objects
- caches are local to each process
- the bridge is not designed to hold multiple active mechanisms/phases at once
  inside the same process
- the bridge should be treated as rank-local and not thread-parallel unless
  additional synchronization is added
```

This is consistent with the current MPI model, where every flow rank owns the
fields it updates and calls Cantera for local/owned data as needed.

---

## 3. Public C ABI entry points

The bridge exposes C ABI routines for Fortran.

### 3.1 Initialization and species metadata

```text
cantera_init_c
cantera_get_species_count_c
cantera_get_species_name_c
```

Responsibilities:

```text
- trim Fortran/C strings
- load the mechanism
- optionally select a named phase
- create the Cantera transport manager
- expose species metadata to Fortran
```

### 3.2 Transport update

```text
cantera_update_transport_c
```

Inputs:

```text
T, P, Y, solver species names
```

Outputs:

```text
mu
D_k
```

The current FV species operator consumes scalar per-species diffusivities.  The
bridge therefore returns scalar `D_k` values using Cantera's mixture-diffusion
coefficient interface.  This is a natural match for mixture-averaged transport.
It is not yet a full multicomponent Stefan-Maxwell flux coupling.

### 3.3 Thermo update from temperature

```text
cantera_update_thermo_c
```

Inputs:

```text
T, P, Y, T_ref, solver species names
```

Outputs:

```text
h_sens
cp
lambda
rho_thermo
```

This path is useful for initialization and pointwise thermo updates from a
known temperature field.

### 3.4 Temperature recovery from enthalpy

```text
cantera_recover_temperature_from_h_c
```

Inputs:

```text
h_sens, P, Y, T_ref, solver species names
```

Output:

```text
T
```

This performs the HP inversion needed because the solver transports enthalpy,
not temperature.

### 3.5 Combined thermo synchronization

```text
cantera_recover_temperature_and_update_thermo_c
```

Inputs:

```text
h_sens, P, Y, T_ref, solver species names
```

Outputs:

```text
T
cp
lambda
rho_thermo
```

This is the preferred energy-step path when the solver needs both temperature
recovery and refreshed thermodynamic properties.  It avoids a separate HP
temperature recovery pass followed by a second TPY property-update pass.

### 3.6 Species sensible enthalpies

```text
cantera_species_sensible_enthalpies_c
```

Outputs species sensible enthalpies:

```text
h_k(T,Y,p0) - h_k(T_ref,Y,p0)
```

These are used by the optional species-enthalpy diffusion correction:

```text
-div(sum_k h_k J_k)
```

### 3.7 Cache diagnostics

```text
cantera_get_cache_stats_c
```

Returns bridge-local cache statistics for:

```text
- transport update
- combined thermo sync
- bulk species-enthalpy calls
- point species-enthalpy calls
```

The counters are stored as floating-point values so the Fortran/MPI side can
reduce them with the same double-precision reduction path used by other
diagnostics.

---

## 4. Named Cantera phase loading

The bridge initialization accepts:

```text
cantera_mech_file
cantera_phase_name
```

Blank phase name:

```text
Cantera::newSolution(mech_file)
```

Nonblank phase name:

```text
Cantera::newSolution(mech_file, phase_name)
```

This lets a case choose the thermodynamic and transport model by selecting a
phase from a YAML file.

Examples of valid phase-level choices include:

```yaml
thermo: ideal-gas
transport: mixture-averaged
```

```yaml
thermo: Peng-Robinson
transport: high-pressure
```

```yaml
thermo: Redlich-Kwong
transport: high-pressure-Chung
```

The solver-level setting:

```fortran
density_eos = "cantera"
```

means:

```text
use density from the selected Cantera phase
```

It does not mean a specific EOS.  The actual EOS is the `thermo:` model in the
selected Cantera phase.

---

## 5. Thermodynamic pressure and state inputs

The bridge receives pressure through arrays passed from Fortran.  In normal
solver usage these pressures represent the thermodynamic/background pressure
used for Cantera state calls:

```text
p0 = params%background_press
```

This pressure is distinct from the hydrodynamic projection pressure field used
by the low-Mach pressure solve.

The bridge assumes that the Fortran side passes the intended thermodynamic
pressure for the property evaluation.  High-pressure EOS cases should therefore
set the solver's thermodynamic pressure consistently with the chosen Cantera
phase and mechanism.

---

## 6. Sensible enthalpy convention

The solver uses sensible enthalpy relative to a reference temperature:

```text
h_sens(T,Y,p0) = h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)
```

Temperature recovery reverses this operation:

```text
h_abs_target = h_sens + h_abs(T_ref,Y,p0)
setState_HP(h_abs_target, p0)
```

This convention matters because non-reacting species mixing should not create
artificial heat release from formation enthalpy.  The transported state is the
sensible mixture enthalpy, while Cantera's absolute enthalpy is used internally
only to perform consistent state updates.

The reference temperature argument `T_ref` must remain consistent across:

```text
cantera_update_thermo_c
cantera_recover_temperature_from_h_c
cantera_recover_temperature_and_update_thermo_c
cantera_species_sensible_enthalpies_c
```

---

## 7. Composition handling

The Fortran solver may transport only a subset of the species in the selected
Cantera phase.  The bridge therefore maps solver species names to Cantera
indices and builds a full Cantera mass-fraction vector for each point.

Current rule:

```text
- solver-provided species mass fractions are inserted by name
- negative transported mass fractions are clipped to zero before Cantera calls
- remaining mass fraction is assigned to N2 when N2 exists
- if no valid mass is present, fall back to N2 if available or the first species
```

This rule supports:

```text
- passive species subsets
- no-species thermo paths
- bath/default composition for non-reacting runs
```

However, for fixed-composition boundary states, the Fortran side should pass the
boundary composition when evaluating boundary enthalpy or boundary density.  The
bridge can evaluate the property once the correct `Y` vector is supplied.

---

## 8. Density semantics

The bridge computes thermodynamic density as:

```text
rho_thermo = gas->density()
```

How that density is used is decided by the Fortran solver mode.

Constant-density mode:

```text
transport%rho = params%rho
energy%rho_thermo = diagnostic / output / validation property
```

Variable-density low-Mach mode:

```text
transport%rho <- energy%rho_thermo
```

Thus, the bridge is the source of Cantera thermodynamic density, but it does not
itself decide whether that density is active in the flow solve.

---

## 9. Transport-model interpretation

The selected Cantera phase owns the transport model.  The bridge creates the
transport manager through:

```cpp
trans = sol->transport();
```

The current bridge exposes scalar transport quantities to Fortran:

```text
mu
lambda
D_k
```

Important implication:

```text
A Cantera phase may be configured as multicomponent or high-pressure transport,
but the current FV species operator still receives scalar D_k values.  A full
multicomponent finite-volume diffusion operator would require additional bridge
outputs and Fortran-side flux assembly.
```

For current validation of species diffusion, mixture-averaged transport is the
closest conceptual match to the solver's scalar-diffusivity species operator.

High-pressure or Chung transport phases are useful for EOS/property testing when
the selected phase supports all queried properties.

---

## 10. Caching strategy

The bridge includes conservative process-local caches to avoid repeated Cantera
work when state variables have not changed.

Current cached paths include:

```text
- transport properties
- combined thermo sync
- species sensible enthalpies
```

Cache dependency rules:

```text
transport cache depends on:
  T, P, Y

combined thermo-sync cache depends on:
  h_sens, P, Y

species-enthalpy cache depends on:
  T, P, Y
```

A cache hit may reuse dependent properties:

```text
T
cp
lambda
rho_thermo
mu
D_k
h_k sensible
```

A cache hit must not reinterpret or overwrite the transported solver state.  The
Fortran arrays remain the source of truth.

---

## 11. Error behavior

The bridge currently treats Cantera failures as fatal:

```text
- print Cantera error or standard exception
- exit(1)
```

This is appropriate for early validation and batch runs where a bad mechanism,
bad phase, invalid state, or unsupported transport model should stop the run
rather than silently continue.

Future hardening could route errors through a Fortran/MPI fatal-error wrapper so
all ranks shut down cleanly with a consistent diagnostic message.

---

## 12. Current limitations

The current bridge is intentionally compact and does not yet provide:

```text
- reaction source terms
- heat-release source terms
- full multicomponent diffusion flux matrices
- Soret/Dufour transport coupling
- explicit boundary-state helper APIs
- multiple simultaneous Cantera phases active in one process
- nonfatal recovery from Cantera HP failures
- thread-safe independent bridge contexts
```

Some of these may never be needed, but they should remain explicit when planning
future chemistry, high-pressure, or multicomponent transport work.

---

## 13. Future development assessment

### 13.1 Near term

```text
- keep this document synchronized with src/cantera_interface.cpp
- add clearer startup reporting of selected thermo and transport models
- expose or log the selected phase's thermo_model and transport_model if practical
- add mechanism/phase validation utility tests
- route bridge fatal errors through MPI-aware solver error handling
```

### 13.2 Short-to-medium term

```text
- add boundary-state helper support:
    rho_b = rho(T_b, Y_b, p0)
    h_b   = h(T_b, Y_b, p0)

- use that support for density-aware variable-density inlet/outlet logic

- add mass-flow boundary support on the Fortran side using boundary EOS density

- extend diagnostics for Cantera cache hit/miss behavior by call site
```

### 13.3 Medium-to-long term

```text
- add reaction-rate source terms and heat release after non-reacting
  variable-density validation is robust

- add a true multicomponent diffusion interface if the FV species operator is
  upgraded from scalar D_k diffusion to coupled Stefan-Maxwell fluxes

- consider a bridge-context object model if future workflows require multiple
  mechanisms/phases active in the same process

- consider optional non-Cantera ideal-gas EOS only if it provides a clear
  performance or dependency benefit
```

---

## 14. Status summary

```text
Mechanism and phase loading:
  active, including optional named phase selection

Cantera thermo sync:
  active and central to energy transport

Cantera density:
  diagnostic in constant-density mode
  active density source in guarded variable-density mode

Cantera transport:
  active for scalar properties consumed by current FV operators

Chemistry:
  not yet active

Full multicomponent FV diffusion:
  not yet implemented

Boundary EOS support:
  needed for production-quality variable-density inlet/outlet and mass-flow BCs
```