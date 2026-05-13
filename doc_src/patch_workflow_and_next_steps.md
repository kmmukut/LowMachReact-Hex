title: Patch Workflow and Next Steps

# Patch Workflow and Next Steps

## Patch workflow

Continue using one small Python patch per logical change.

Each patch should:

```text
- support --dry-run and --apply
- create timestamped backups under .backups/YYYYMMDD_HHMMSS/
- fail loudly when anchors are missing or ambiguous
- be idempotent where practical
- print a clear summary of changed/skipped files
- list build and validation next steps
```

## Recommended immediate sequence

1. Keep Patch 019 profiling updates.
2. Keep Patch 020 combined Cantera thermo-sync optimization.
3. Inspect `git diff`.
4. Build release and debug.
5. Run energy-disabled, constant-property energy, Cantera thermo without species, and Cantera thermo with species cases.
6. Compare profiling output before/after thermo-sync optimization.
7. Only then proceed to a radiation coupling stub, manufactured-qrad patch, or lower-risk diagnostics such as Cantera cache hit/miss counters.

## Current design decision

The solver uses Option A.

When species transport changes `Y` before energy transport, preserve transported enthalpy:

```text
h_after_species = h_before_species
T_after_species = T(h_after_species, Y_new, p0)
```

Do not reset `h` from the old temperature after species changes:

```text
do not use: h = h(T_old, Y_new, p0)
```

The current Cantera thermo-sync optimization preserves this convention. It may combine temperature recovery and property refresh into one Cantera pass, but it must not alter the transported `h` field.
