title: Profiling Guide

# Profiling Guide

This guide explains how to enable and interpret the LowMachReact-Hex terminal profiler.

The profiler is intended to answer practical performance questions:

```text
Where did the runtime go?
Did a patch actually make the intended region faster?
Is the bottleneck projection, species, energy, Cantera, MPI, diagnostics, or output?
```

The current profiler is MPI-aware and can print both a flat timing table and a nested call tree.

## Enabling profiling

Enable profiling in `case.nml`:

```fortran
&profiling_input
  enable_profiling = .true.
  nested_profiling = .true.
/
```

Recommended usage:

```text
enable_profiling = .true.
nested_profiling = .true.
```

Use nested profiling for development and optimization. The nested tree is the safest way to interpret parent/child costs.

For production timing runs, use the same:

```text
- case
- mesh
- number of MPI ranks
- compiler/build mode
- nsteps
- output_interval
- profiling settings
```

when comparing before/after results.

## Report structure

A typical profiling report contains two sections:

```text
PERFORMANCE PROFILING REPORT
NESTED PROFILING REPORT
```

The flat report lists all recorded timers. The nested report shows the parent/child region tree.

The report prints average timing across MPI ranks, plus minimum and maximum rank times for each timer.

Important rule:

```text
Timings are inclusive.
Flat rows are not additive when nested timers are enabled.
```

For example, if the report shows:

```text
Energy_Transport
  Energy_Cantera_PreSync
  Energy_Flux_Update
  Energy_Cantera_PostSync
```

then the child timer costs are already included inside `Energy_Transport`.

Do not add `Energy_Transport` and its child timers together.

## Main top-level regions

Current top-level regions include:

```text
Total_Simulation
Transport_Update
Projection_Step
Species_Transport
Energy_Transport
CFL_Update
Flow_Diagnostics
Diagnostics_Write_Flow
Diagnostics_Write_Energy
Output_Write_VTU
```

### `Total_Simulation`

Overall wall time for the full simulation.

Use this as the denominator when comparing total speedup.

### `Transport_Update`

Updates transport properties on `transport_update_interval`.

This may include:

```text
Transport_Setup
Transport_Cantera_Call
Transport_Unpack
Transport_Exchange
Transport_Exchange_Diff
```

This region controls Cantera transport properties such as mixture viscosity and species diffusivity, not the energy-side Cantera thermo sync.

### `Projection_Step`

Advances momentum and the pressure projection.

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

If `Projection_Step` dominates, inspect `Projection_PCG` first.

### `Species_Transport`

Advances passive species transport.

If this dominates, inspect species advection/diffusion work, species count, and species halo exchange.

### `Energy_Transport`

Advances passive sensible enthalpy transport and thermodynamic synchronization.

For Cantera thermo runs after the combined thermo-sync optimization, common children are:

```text
Energy_Exchange_H
Energy_Cantera_PreSync
Energy_PreFlux_Exchange
Energy_Flux_Update
Energy_Cantera_PostSync
Energy_Final_Exchange
```

For constant-cp energy runs, the Cantera sync timers should be absent.

### `Diagnostics_Write_Flow`

Writes flow diagnostics.

This should usually be small.

### `Diagnostics_Write_Energy`

Writes energy diagnostics.

This should usually be small.

### `Output_Write_VTU`

Writes VTU/PVTU output.

If this dominates, reduce output frequency or inspect output array size.

## Energy profiling interpretation

The energy region separates thermodynamic synchronization from the finite-volume energy update.

### `Energy_Cantera_PreSync`

Synchronizes dependent thermo state before the energy flux update:

```text
(T, cp, lambda, rho_thermo) = sync(h, Y, p0)
```

This pre-flux sync is needed when transported species may have changed composition before the energy update.

For species-disabled Cantera thermo runs, this timer may be absent because composition has not changed and the previous post-flux sync remains valid.

### `Energy_PreFlux_Exchange`

Synchronizes off-rank temperature and thermal-conductivity values needed by the conductive flux stencil:

```text
T
lambda
```

If this timer is large, the cost is usually halo exchange.

### `Energy_Flux_Update`

The finite-volume enthalpy update itself.

This includes:

```text
upwind h advection
Fourier conduction driven by grad(T)
qrad/rho source
```

If this timer is small, the energy stencil is not the bottleneck.

### `Energy_Cantera_PostSync`

Synchronizes dependent thermo state after the enthalpy update:

```text
(T, cp, lambda, rho_thermo) = sync(h_new, Y, p0)
```

This prepares output and the next step.

### `Energy_Final_Exchange`

Synchronizes updated energy fields for output and the next step.

Typical fields include:

```text
h
T
lambda
```

## Cantera thermo-sync optimization

The current energy-side Cantera thermo path uses a combined sync operation where possible:

```text
(T, cp, lambda, rho_thermo) = sync(h, Y, p0)
```

This is logically equivalent to:

```text
recover T from h,Y,p0
refresh cp/lambda/rho_thermo from T,Y,p0
```

but avoids a redundant second full Cantera pass.

The sync operation must preserve transported enthalpy:

```text
h_after_sync = h_before_sync
```

The combined thermo-sync path may use a conservative cache. The cache key is the transported thermodynamic state:

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

## MPI profiling interpretation

`MPI_Communication` may appear both as a flat timer and as a child of other regions.

Common locations:

```text
Projection_PCG -> MPI_Communication
Transport_Exchange -> MPI_Communication
Species_Transport -> MPI_Communication
Energy_PreFlux_Exchange -> MPI_Communication
Energy_Final_Exchange -> MPI_Communication
```

### Pressure-solver MPI

If the nested report shows:

```text
Projection_PCG
  MPI_Communication
```

then pressure-solver communication is a major cost.

This may come from:

```text
halo exchanges in pressure matvecs
global dot-product reductions
pressure correction synchronization
```

Pressure/PCG MPI optimizations can easily affect solver correctness, so treat them as high-risk changes.

### Energy MPI

If the nested report shows:

```text
Energy_PreFlux_Exchange
  MPI_Communication
Energy_Final_Exchange
  MPI_Communication
```

then the energy path is paying for halo exchanges of temperature, conductivity, enthalpy, or related fields.

Possible future optimizations include packed multi-scalar halo exchange, but only after correctness tests pass.

### MPI time can indicate load imbalance

High MPI time does not always mean the communication routine itself is inefficient.

It can also mean some ranks are waiting for slower ranks.

Possible causes:

```text
uneven owned-cell count
uneven Cantera cache misses
boundary-heavy partitions
composition-gradient-heavy partitions
node noise
```

Useful diagnostics for future optimization:

```text
owned cells per rank
Cantera thermo cache hits/misses per rank
pressure iterations per step
MPI min/max timing spread
```

## How to compare two runs

Use this workflow when evaluating a patch:

1. Run the same case before and after the change.
2. Use the same number of MPI ranks.
3. Use the same build type.
4. Use the same `nsteps`.
5. Use the same `output_interval`.
6. Save the full terminal output from both runs.
7. Compare `Total_Simulation`.
8. Compare top-level regions.
9. Compare nested child timers inside the changed region.
10. Confirm diagnostics remain physically equivalent.

Do not judge a performance patch only from total runtime. Always check that the intended region changed.

Example:

```text
If a Cantera thermo patch was applied:
    compare Energy_Cantera_PreSync/PostSync or old recover/refresh timers.

If an output patch was applied:
    compare Output_Write_VTU.

If a pressure-solver patch was applied:
    compare Projection_PCG, pressure iterations, and divergence diagnostics.
```

## Common profiling patterns

### Cantera thermo dominates

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
add cache hit/miss diagnostics
```

### Energy stencil dominates

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
```

### Projection dominates

Typical tree:

```text
Projection_Step
  Projection_PCG
    MPI_Communication
    Pressure_Matvec
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
inspect MPI reductions
do not change PCG communication casually
```

### Output dominates

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
reduce number of output arrays
inspect filesystem performance
```

## Common mistakes

### Mistake 1: Adding flat rows together

Do not add parent and child timers together.

Wrong:

```text
Energy_Transport + Energy_Cantera_PreSync + Energy_Cantera_PostSync
```

Correct:

```text
Energy_Transport already includes its children.
```

### Mistake 2: Comparing different cases

Changing mesh size, rank count, timestep count, or output interval invalidates a direct performance comparison.

### Mistake 3: Optimizing tiny timers first

Focus on large parent regions first.

For example, if:

```text
Projection_PCG = 30%
Output_Write_VTU = 0.5%
```

then output is not the immediate bottleneck.

### Mistake 4: Treating MPI time as automatically bad

MPI time may reflect communication overhead, synchronization cost, or load imbalance.

Check min/max timing spread before assuming the communication routine itself is the problem.

### Mistake 5: Ignoring correctness after speedup

Every performance change must preserve:

```text
max/rms divergence
net boundary flux
temperature bounds
enthalpy diagnostics
species boundedness
sum_Y behavior
```

Speedup without physical equivalence is not a successful optimization.

## Recommended profiling checkpoints

For the current staged solver, useful profiling checkpoints are:

```text
1. Energy disabled baseline
2. Constant-property energy
3. Cantera thermo without transported species
4. Cantera thermo with transported species
5. Cantera species diffusivity with fixed-Re flow
6. Non-reacting counterflow development case
```

Save representative reports for each stage. They make future regressions much easier to identify.

## Practical acceptance criteria for profiling patches

A profiling or performance patch should satisfy:

```text
- diagnostics remain equivalent
- output fields remain equivalent
- Total_Simulation improves or the target timer improves
- new timers make interpretation clearer
- profiler overhead remains small compared with total runtime
```

For documentation or timer-name-only patches, numerical results should be exactly unchanged apart from normal floating-point and output-order effects.
