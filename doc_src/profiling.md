title: Profiling Guide


# Profiling Guide

This guide explains how to enable and interpret the LowMachReact-Hex terminal
profiler and Cantera cache statistics.

The profiler is intended to answer practical performance questions:

```text
Where did the runtime go?
Did a patch actually make the intended region faster?
Is the bottleneck projection, species, energy, Cantera, MPI, diagnostics, or output?
Did the speedup preserve the numerical diagnostics?
```

The profiler is MPI-aware and can print both a flat timing table and a nested
call tree.

---

## 1. Enabling profiling

Enable profiling in `case.nml`:

```fortran
&profiling_input
  enable_profiling = .true.
  nested_profiling = .true.
  write_cantera_cache_stats = .false.
  variable_density_debug = .false.
/
```

Recommended development setting:

```fortran
enable_profiling = .true.
nested_profiling = .true.
write_cantera_cache_stats = .false.
variable_density_debug = .false.
```

Use `nested_profiling = .true.` for development and optimization.  The nested
tree is the safest way to interpret parent/child costs.

Use `write_cantera_cache_stats = .true.` only when you need Cantera cache
diagnostics.

Use `variable_density_debug = .true.` only when debugging variable-density
stability or state issues.  It is not a normal profiling option and can add
extra terminal output.

---

## 2. What to keep fixed when comparing runs

For a fair before/after performance comparison, keep these identical:

```text
case
mesh
number of MPI ranks
compiler and build mode
nsteps
dt or dynamic timestep settings
output_interval
physics flags
Cantera mechanism and phase
profiling settings
diagnostics/output settings
machine allocation where possible
```

Changing output cadence, mesh size, number of species, Cantera phase, MPI rank
count, or build type invalidates a direct timing comparison.

For noisy systems, run more than once and compare trends, not one isolated
number.

---

## 3. Report structure

A typical profiling report contains:

```text
PERFORMANCE PROFILING REPORT
NESTED PROFILING REPORT
```

The flat report lists recorded timers.  The nested report shows parent/child
regions.

The report prints timing aggregated across MPI ranks, including average, minimum,
and maximum rank timing where supported.

Important rule:

```text
Timings are inclusive.
Flat rows are not additive when nested timers are enabled.
```

Example:

```text
Energy_Transport
  Energy_Cantera_PreSync
  Energy_Flux_Update
  Energy_Cantera_PostSync
```

The child timers are already included inside `Energy_Transport`.  Do not add
`Energy_Transport` and its children together.

Use the nested report to understand where time inside a parent region went.

---

## 4. Main top-level regions

Common top-level regions include:

```text
Total_Simulation
CFL_Update
Transport_Update
Projection_Step
Species_Transport
Energy_Transport
Flow_Diagnostics
Diagnostics_Write_Flow
Diagnostics_Write_Energy
Output_Write_VTU
```

Not every timer appears in every run.  Timers are gated by enabled physics,
diagnostics, and output settings.

### 4.1 `Total_Simulation`

Overall wall time for the simulated run.

Use this as the denominator for total speedup.

### 4.2 `CFL_Update`

Computes CFL and dynamic timestep information when active.

This should usually be small unless the case is very large or reductions are
costly.

### 4.3 `Transport_Update`

Refreshes transport properties on `transport_update_interval`.

This may include:

```text
Transport_Setup
Transport_Cantera_Call
Transport_Unpack
Transport_Exchange
Transport_Exchange_Diff
```

This region covers Cantera-backed transport properties such as:

```text
mu
D_k
```

It is separate from the energy-side Cantera thermo sync.

### 4.4 `Projection_Step`

Advances momentum and pressure projection.

Common children include:

```text
Projection_Momentum_RHS
Projection_AB2
Projection_Predict_Flux
Projection_Poisson_RHS
Projection_PCG
Projection_Pressure_Update
Projection_Correction
Projection_Diagnostics
```

If `Projection_Step` dominates, inspect:

```text
Projection_PCG
pressure iteration count
pressure residual
MPI_Communication under Projection_PCG
```

Projection optimizations are high-risk because they can affect continuity,
boundary compatibility, and pressure nullspace handling.

### 4.5 `Species_Transport`

Advances species transport when enabled.

Costs usually scale with:

```text
number of transported species
number of faces/cells
diffusion model
halo exchange cost
boundary condition complexity
```

In variable-density mode, species transport uses the conservative `rho*Y`
branch, so performance can differ from constant-density species transport.

### 4.6 `Energy_Transport`

Advances sensible enthalpy and thermo synchronization when enabled.

Common Cantera-thermo children include:

```text
Energy_Exchange_H
Energy_Cantera_PreSync
Energy_PreFlux_Exchange
Energy_Flux_Update
Energy_Cantera_PostSync
Energy_Final_Exchange
Energy_Cantera_SpeciesH
```

Not all children appear in all modes.

For constant-cp energy runs, Cantera sync timers should be absent.

For species-disabled Cantera thermo runs, `Energy_Cantera_PreSync` may be absent
or much smaller because composition has not changed and the previous post-sync
state may already be valid.

### 4.7 `Diagnostics_Write_Flow`

Writes flow diagnostics.

Usually small.  If it is large, check output cadence, filesystem performance,
and whether very detailed diagnostics are enabled.

### 4.8 `Diagnostics_Write_Energy`

Writes energy diagnostics.

Usually small, but variable-density energy validation may add more reductions
and CSV columns.

### 4.9 `Output_Write_VTU`

Writes VTU/PVTU/PVD visualization files.

If this dominates:

```text
increase output_interval
reduce number of output arrays if appropriate
check filesystem performance
avoid unnecessary debug fields in production runs
```

---

## 5. Energy profiling interpretation

The energy region separates thermodynamic synchronization from the finite-volume
enthalpy update.

### 5.1 `Energy_Cantera_PreSync`

Synchronizes dependent thermo state before the energy flux update:

```text
(T, cp, lambda, rho_thermo) = sync(h, Y, p0)
```

This pre-flux sync is needed when transported species may have changed
composition before the energy update.

For species-disabled Cantera thermo runs, this timer may be absent because the
previous post-flux sync remains valid.

### 5.2 `Energy_PreFlux_Exchange`

Synchronizes fields needed by the conductive/energy flux stencil, commonly:

```text
T
lambda
```

If this timer is large, the cost is usually halo exchange or waiting for slow
ranks.

### 5.3 `Energy_Flux_Update`

Finite-volume enthalpy update.

This includes the enabled subset of:

```text
upwind h advection
Fourier conduction driven by grad(T)
qrad source contribution
species-enthalpy diffusion contribution
rho*h conservative update in variable-density mode
```

If this dominates, inspect face-loop cost, boundary work, memory access pattern,
and species-enthalpy diffusion work.

### 5.4 `Energy_Cantera_PostSync`

Synchronizes dependent thermo state after the enthalpy update:

```text
(T, cp, lambda, rho_thermo) = sync(h_new, Y, p0)
```

This prepares the next step, output, diagnostics, and variable-density density
sync.

### 5.5 `Energy_Final_Exchange`

Synchronizes updated energy fields for output and later operators.

Typical fields include:

```text
h
T
lambda
```

### 5.6 `Energy_Cantera_SpeciesH`

Appears when species-enthalpy diffusion is enabled.

This timer covers Cantera bridge work used to compute species sensible enthalpies
for the correction:

```text
-div(sum_k h_k J_k)
```

The bridge uses a conservative cache keyed on:

```text
T
p0
transported species mass fractions
```

If this timer dominates, compare cache hit/miss statistics and consider whether
species-enthalpy diffusion is necessary for the current validation run.

---

## 6. Cantera thermo-sync optimization

The energy-side Cantera thermo path uses a combined sync where possible:

```text
(T, cp, lambda, rho_thermo) = sync(h, Y, p0)
```

This is logically equivalent to:

```text
1. recover T from h,Y,p0
2. refresh cp, lambda, rho_thermo from T,Y,p0
```

but avoids a redundant second full Cantera pass.

The sync must preserve transported enthalpy:

```text
h_after_sync = h_before_sync
```

The cache key for combined thermo sync is the transported thermodynamic state:

```text
h_sensible
Y
p0
```

The cache may reuse only dependent thermodynamic fields:

```text
T
cp
lambda
rho_thermo
```

It must not overwrite or reinterpret the transported `h` field.

---

## 7. Cantera transport versus Cantera thermo

Do not confuse the two Cantera-heavy paths.

### 7.1 Transport Cantera path

Appears under:

```text
Transport_Update
```

Provides:

```text
mu
D_k
```

Controlled by:

```text
fluid_input enable_cantera
species_input enable_cantera
transport_update_interval
```

### 7.2 Energy thermo Cantera path

Appears under:

```text
Energy_Transport
```

Provides:

```text
T(h,Y,p0)
cp
lambda
rho_thermo
species sensible enthalpies when needed
```

Controlled by:

```text
enable_energy
enable_cantera_thermo
enable_species_enthalpy_diffusion
```

### 7.3 Variable-density implication

In constant-density mode:

```text
rho_thermo is diagnostic
```

In variable-density mode:

```text
transport%rho <- energy%rho_thermo
```

So an energy thermo-sync performance change can affect variable-density runs more
directly than constant-density runs.

---

## 8. Cantera cache statistics

Set:

```fortran
write_cantera_cache_stats = .true.
```

to print the final `CANTERA CACHE STATISTICS` block.  The counters are summed
across flow ranks.

The report may include:

```text
Transport_mu_Dk
Energy_ThermoSync
SpeciesH_Bulk
SpeciesH_Point
```

Meaning:

| Counter group | Meaning |
|---|---|
| `Transport_mu_Dk` | Cantera transport calls for viscosity and species diffusivities. |
| `Energy_ThermoSync` | Combined energy thermo-sync calls. |
| `SpeciesH_Bulk` | Bulk cell species sensible-enthalpy calls for species-enthalpy diffusion. |
| `SpeciesH_Point` | Point or boundary species sensible-enthalpy calls. |

Typical reported values include:

```text
calls
states processed
cache hits
cache misses
hit percentage
```

Interpretation:

```text
high hit rate + high Cantera timer:
  each miss may still be expensive or many states are processed

low hit rate:
  state is changing frequently, cache tolerance/keying may not help, or the case
  has strong spatial variation

high SpeciesH_Point:
  boundary/point species-enthalpy evaluations may be significant
```

Use cache statistics as performance diagnostics, not physics validation metrics.

---

## 9. MPI profiling interpretation

`MPI_Communication` may appear as a flat timer and as a child of other regions.

Common locations:

```text
Projection_PCG -> MPI_Communication
Transport_Exchange -> MPI_Communication
Species_Transport -> MPI_Communication
Energy_PreFlux_Exchange -> MPI_Communication
Energy_Final_Exchange -> MPI_Communication
Output_Write_VTU -> MPI_Communication
```

### 9.1 Pressure-solver MPI

If the nested report shows:

```text
Projection_PCG
  MPI_Communication
```

then pressure-solver communication is a major cost.

Possible sources:

```text
pressure matvec halo exchange
global dot-product reductions
pressure correction synchronization
load imbalance visible as waiting time
```

Pressure/PCG MPI optimizations can easily affect solver correctness.  Treat them
as high-risk numerical/infrastructure changes.

### 9.2 Energy MPI

If the nested report shows:

```text
Energy_PreFlux_Exchange
  MPI_Communication
Energy_Final_Exchange
  MPI_Communication
```

then the energy path is paying for halo exchanges of temperature, conductivity,
enthalpy, or related fields.

Possible future optimization:

```text
packed multi-scalar halo exchange
```

but only after correctness tests pass.

### 9.3 MPI time can indicate load imbalance

High MPI time does not automatically mean the communication routine itself is
inefficient.

It can also mean ranks are waiting for slower ranks.

Possible causes:

```text
uneven owned-cell count
uneven Cantera cache misses
boundary-heavy partitions
composition-gradient-heavy partitions
filesystem contention
node noise
```

Useful supporting diagnostics:

```text
owned cells per rank
Cantera cache hits/misses per rank if available
pressure iterations per step
MPI min/max timing spread
timer max/avg ratio
```

---

## 10. How to compare two runs

Use this workflow:

```text
1. Run the baseline case.
2. Save terminal output and diagnostics.
3. Apply the patch.
4. Run the same case with the same settings.
5. Compare Total_Simulation.
6. Compare top-level regions.
7. Compare nested child timers inside the changed region.
8. Compare Cantera cache statistics if Cantera was affected.
9. Confirm diagnostics remain physically equivalent.
```

Examples:

```text
Cantera thermo patch:
  compare Energy_Cantera_PreSync, Energy_Cantera_PostSync, Energy_ThermoSync cache stats

species-enthalpy diffusion patch:
  compare Energy_Cantera_SpeciesH and SpeciesH_Bulk/Point stats

output patch:
  compare Output_Write_VTU and check PVTU schema

pressure-solver patch:
  compare Projection_PCG, pressure iterations, projection residuals, and continuity diagnostics

variable-density source/projection patch:
  compare divu_minus_S_projection_*, conservative continuity metrics, and Total_Simulation
```

Do not judge a performance patch only from total runtime.  Always verify that
the intended timer changed and that numerical diagnostics remain acceptable.

---

## 11. Common profiling patterns

### 11.1 Cantera thermo dominates

Typical tree:

```text
Energy_Transport
  Energy_Cantera_PreSync
  Energy_Cantera_PostSync
```

Interpretation:

```text
The bottleneck is thermodynamic recovery/property refresh.
```

Possible responses:

```text
check combined thermo-sync path
check cache behavior
avoid redundant pre-sync when species are disabled
profile mechanism/phase complexity
```

### 11.2 Species sensible enthalpy dominates

Typical tree:

```text
Energy_Transport
  Energy_Cantera_SpeciesH
```

Interpretation:

```text
Species-enthalpy diffusion property evaluation is expensive.
```

Possible responses:

```text
inspect SpeciesH_Bulk and SpeciesH_Point cache stats
verify whether enable_species_enthalpy_diffusion is needed for this run
look for repeated boundary/point evaluations
```

### 11.3 Energy stencil dominates

Typical tree:

```text
Energy_Transport
  Energy_Flux_Update
```

Interpretation:

```text
The finite-volume enthalpy stencil is the bottleneck.
```

Possible responses:

```text
inspect face-loop work
inspect boundary enthalpy evaluations
check memory access pattern
check whether variable-density rho*h branch is active
```

### 11.4 Projection dominates

Typical tree:

```text
Projection_Step
  Projection_PCG
    MPI_Communication
```

Interpretation:

```text
The pressure solve is the bottleneck.
```

Possible responses:

```text
inspect pressure iterations
inspect pressure tolerance
inspect preconditioner
inspect matrix coefficient behavior in variable-density mode
do not change PCG communication casually
```

### 11.5 Output dominates

Typical tree:

```text
Output_Write_VTU
```

Interpretation:

```text
Visualization output is expensive.
```

Possible responses:

```text
increase output_interval
reduce optional/debug output arrays when appropriate
inspect filesystem performance
remove stale output before schema-changing tests
```

### 11.6 Diagnostics dominate

Typical tree:

```text
Diagnostics_Write_Flow
Diagnostics_Write_Energy
```

Interpretation:

```text
Diagnostics are expensive or too frequent.
```

Possible responses:

```text
increase output_interval
disable verbose debug diagnostics
avoid per-cell audit outputs except during targeted debugging
```

---

## 12. Common mistakes

### Mistake 1: Adding flat rows together

Wrong:

```text
Energy_Transport + Energy_Cantera_PreSync + Energy_Cantera_PostSync
```

Correct:

```text
Energy_Transport already includes its children.
```

### Mistake 2: Comparing different cases

Changing mesh size, rank count, timestep count, build type, output interval, or
Cantera mechanism invalidates a direct comparison.

### Mistake 3: Optimizing tiny timers first

Focus on large parent regions.

Example:

```text
Projection_PCG = 30%
Output_Write_VTU = 0.5%
```

In that case, output is not the immediate bottleneck.

### Mistake 4: Treating MPI time as automatically bad

MPI time may reflect communication overhead, synchronization cost, load
imbalance, filesystem contention, or waiting for slow ranks.

Check min/max timing spread before assuming the communication routine itself is
the problem.

### Mistake 5: Ignoring correctness after speedup

Every performance change must preserve relevant diagnostics:

```text
constant-density:
  max/rms divergence
  net boundary flux
  pressure iterations/residuals

species:
  boundedness
  sum_Y
  species integrals

energy:
  temperature bounds
  enthalpy diagnostics
  energy closure metrics

variable-density:
  divu_minus_S_projection_*
  conservative continuity metrics
  rho*h closure/reconciliation metrics
```

Speedup without physical equivalence is not a successful optimization.

### Mistake 6: Profiling with excessive output

Frequent VTU writes can dominate runtime and hide solver costs.

Use a larger `output_interval` when measuring solver kernels.

### Mistake 7: Leaving `variable_density_debug` on

Verbose variable-density debug output can distort terminal-heavy runs and make
profiling harder to read.  Use it only for targeted debugging.

---

## 13. Recommended profiling checkpoints

Keep representative profiling reports for these stages:

```text
1. Flow-only constant-density baseline
2. Constant-property energy
3. Cantera thermo without transported species
4. Cantera thermo with transported species
5. Cantera species diffusivity with fixed-Re flow
6. Species-enthalpy diffusion enabled
7. Guarded non-reacting variable-density low-Mach case
8. Representative EOS/pressure variants if using real-gas mechanisms
```

These checkpoints make performance regressions easier to identify as the solver
evolves.

---

## 14. Practical acceptance criteria for performance patches

A profiling or performance patch should satisfy:

```text
- debug build passes
- release build passes
- diagnostics remain equivalent for affected modes
- output fields remain equivalent unless the patch intentionally changes output
- Total_Simulation improves or the targeted timer improves
- new timers make interpretation clearer
- profiler overhead remains small relative to total runtime
- documentation is updated if timer names or output behavior change
```

For documentation-only or timer-name-only patches, numerical results should be
unchanged apart from normal floating-point and output-order effects.

---

## 15. Related diagnostics

Profiling tells where time is spent.  It does not by itself prove numerical
correctness.

Use these diagnostics with profiling:

```text
diagnostics.csv
energy_diagnostics.csv
species_energy_conservation.csv
species_integrals.csv
enthalpy_energy_budget.csv
variable_density_diagnostics.csv
variable_density_compatibility.csv
variable_density_transport_conservation.csv
variable_density_continuity_residual.csv
```

For variable-density validation, primary pass/fail metrics are documented in:

```text
doc_src/validation_metrics.md
```

For output-file locations, see:

```text
doc_src/output_layout.md
```

For field meanings in ParaView, see:

```text
doc_src/paraview_output_fields.md
```

---

## 16. Maintenance rule

Keep this guide focused on profiling interpretation and performance workflow.

Do not append long historical patch notes here.  This guide should answer:

```text
How do I enable profiling?
What does each timer mean?
How do I compare two runs?
How do I interpret Cantera cache statistics?
How do I avoid misleading timing comparisons?
Which correctness diagnostics must remain valid?
```