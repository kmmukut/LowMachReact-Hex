title: Validation Metrics

# Validation Metrics

This document defines the accepted current-state validation metrics for
LowMachReact-Hex. It is intentionally about metric meaning and pass/fail
interpretation. Automation commands live in `validation_automation.md`.

## 1. Scope

The stable baseline remains constant-density projection. The guarded
variable-density validation path is:

```text
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
thermo_update_interval = 1
enable_reactions = .false.
```

`density_eos="cantera"` means active density comes from the selected Cantera
phase. The actual EOS is the YAML phase model.

## 2. Projection

Constant-density target:

$$
\nabla \cdot u = 0 .
$$

Variable-density target:

$$
\nabla \cdot u = S_{projection}.
$$

Primary variable-density projection metrics are the
`divu_minus_S_projection_*` columns. Raw `div(u)` and current-source residuals
are explanatory/source-evolution diagnostics, not the primary projection error.

## 3. Density

Constant-density mode:

```text
transport%rho = params%rho
energy%rho_thermo = diagnostic
```

Variable-density mode:

```text
transport%rho <- energy%rho_thermo
```

Mass diagnostics should use active `transport%rho`.

## 4. Continuity

The conservative continuity target is:

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) = 0 .
$$

Use the conservative residual and coupled-transport audit outputs as primary
continuity checks. Expanded/current-source residuals help diagnose time-level
or source-evolution issues.

## 5. Species

Supported variable-density species transport is non-reacting and conservative
in \(\rho Y_k\). The expected checks are:

```text
0 <= Y_k <= 1
sum_k Y_k ~= 1
species mass changes are consistent with boundary fluxes and diffusion
```

The correction velocity should keep the net diffusive species mass flux near
zero.

## 6. Energy

Energy transports sensible enthalpy `h`, not temperature. Option A is active:
preserve transported `h`; recover \(T = T(h,Y,p_0)\) after species changes.

Primary energy checks are:

| Metric family | Meaning |
|---|---|
| Direct update closure | Checks the discrete enthalpy update actually applied. |
| Operator-consistent budget | Uses the density/time level of the update operator. |
| Output-state reconciliation | Explains differences between output-state `rho*h` and operator-consistent `rho*h`. |

`qrad` is a source hook. Unless a case explicitly fills it, physical radiation
coupling is not being validated.

## 7. Output Sources

Solver-generated diagnostics live under:

```text
<output_dir>/diagnostics/
```

Visualization fields under `<output_dir>/VTK/` are useful for localization, but
CSV diagnostics are the source of truth for pass/fail validation.

## 8. Automation

Use `validation_automation.md` for checker and matrix-runner workflows. Those
tools should implement the metric meanings in this document rather than
defining independent physics criteria.
