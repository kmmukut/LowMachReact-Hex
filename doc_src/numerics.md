title: Numerical Method

# Numerical Method

LowMachReact-Hex currently implements a collocated finite-volume solver for constant-density low-Mach / incompressible flow on hexahedral control volumes. The current staged physics supports optional passive species transport and optional passive sensible-enthalpy transport with Cantera thermodynamic/property evaluation.

The current implementation should be interpreted as a **constant-density projection solver**, not a variable-density reacting low-Mach solver. The density used by the momentum equation, conservative face fluxes, CFL estimate, pressure projection, and diagnostic mass integral is the configured flow density:

$$
\rho_f = \rho_{\text{params}}
$$

Cantera thermodynamic density is stored separately as a diagnostic field:

$$
\rho_{\text{thermo}} = \rho_{\text{Cantera}}(T,Y,p_0)
$$

and is not used in the pressure projection or continuity constraint.

## Enabled and Disabled Physics

The currently supported model contains:

- constant flow / projection density,
- incompressible pressure projection,
- optional passive species mass-fraction transport,
- optional passive sensible-enthalpy transport,
- Fourier heat conduction driven by temperature gradients,
- Cantera evaluation of \(h(T,Y,p_0)\), \(T(h,Y,p_0)\), \(c_p\), \(\lambda\), and diagnostic \(\rho_{\text{thermo}}\),
- radiation-source storage through `qrad`, with `qrad = 0` unless filled by a future coupling.

The following are intentionally not part of the current model:

- variable-density low-Mach divergence constraint,
- reactions,
- reaction heat release,
- species-diffusion enthalpy correction,
- external radiation physics,
- radiation MPI or spectral decomposition.

In particular, the energy equation does **not** yet include:

$$
-\nabla \cdot \left(\sum_k h_k \mathbf{J}_k\right)
$$

where \(h_k\) is the species enthalpy and \(\mathbf{J}_k\) is the diffusive species mass flux. Therefore, composition diffusion is modeled, but enthalpy transport carried by diffusing species is not yet modeled.

## Notation

Cell-centered quantities are stored at control-volume centers. For a cell \(c\):

- \(V_c\) is the cell volume.
- \(f \in \partial c\) denotes a face of cell \(c\).
- \(A_f\) is the face area.
- \(d_f\) is the normal distance used for the face gradient.
- \(\mathbf{n}_{f,c}\) is the unit normal pointing outward from cell \(c\).
- \(F_{f,c}\) is the volumetric face flux oriented outward from cell \(c\).

The stored face flux is oriented from owner cell to neighbor cell. During cell updates, the solver reorients it locally so that positive \(F_{f,c}\) always means flow leaving the current cell.

For a scalar \(\psi\) transported by upwind advection:

$$
\psi_f^{\text{up}} =
\begin{cases}
\psi_c, & F_{f,c} \ge 0, \\
\psi_{\text{other}}, & F_{f,c} < 0.
\end{cases}
$$

Here \(\psi_{\text{other}}\) is either the neighboring-cell value or a boundary value.

## Governing Equations

### Incompressible Flow

The velocity field satisfies the constant-density incompressible constraint:

$$
\nabla \cdot \mathbf{u} = 0
$$

The momentum equation is advanced in fractional-step form:

$$
\frac{\partial \mathbf{u}}{\partial t}
+
\nabla \cdot (\mathbf{u}\mathbf{u})
=
-\frac{1}{\rho_f}\nabla p
+
\nabla \cdot (\nu \nabla \mathbf{u})
+
\mathbf{f}
$$

where \(\rho_f\) is the constant configured flow density, \(\nu\) is the kinematic viscosity used by the flow solver, and \(\mathbf{f}\) is the configured body acceleration.

When Cantera transport is enabled, Cantera may supply mixture viscosity \(\mu\). The flow kinematic viscosity is then computed using the constant projection density,

$$
\nu = \frac{\mu}{\rho_f},
$$

only when variable-viscosity flow is explicitly enabled. Otherwise, validation runs keep the configured constant viscosity.

### Passive Species

For each transported species \(k = 1,\dots,N_s\), the current model advances passive mass fractions:

$$
\frac{\partial Y_k}{\partial t}
+
\nabla \cdot (\mathbf{u}Y_k)
=
\nabla \cdot (D_k \nabla Y_k)
$$

with no reaction source:

$$
\dot{\omega}_k = 0.
$$

The solver applies a diffusive-flux correction so the sum of diffusive species fluxes is zero at each face. This keeps the species set consistent with mass-fraction transport.

### Passive Sensible Enthalpy

The transported thermal variable is mixture sensible enthalpy:

$$
h = h_{\text{sensible}}(T,Y,p_0)
$$

not temperature. The energy equation currently solved is:

$$
\frac{\partial h}{\partial t}
+
\nabla \cdot (\mathbf{u}h)
=
\frac{1}{\rho_f}\nabla \cdot (\lambda \nabla T)
+
\frac{q_{\text{rad}}}{\rho_f}
$$

where \(\lambda\) is thermal conductivity and \(q_{\text{rad}}\) is a volumetric source term in W/m\(^3\).

The conduction operator uses \(\nabla T\), not \(\nabla h\).

## Cantera Thermodynamic Convention

When Cantera thermodynamics are enabled, the solver stores **sensible enthalpy relative to a reference temperature**. For a local composition \(Y\) and pressure \(p_0\):

$$
h(T,Y,p_0)
=
h_{\text{abs}}^{\text{Cantera}}(T,Y,p_0)
-
h_{\text{abs}}^{\text{Cantera}}(T_{\text{ref}},Y,p_0)
$$

where \(T_{\text{ref}}\) is `energy_reference_T`.

Temperature recovery is performed by converting the stored sensible enthalpy back to an absolute Cantera enthalpy target:

$$
h_{\text{target,abs}}
=
h
+
h_{\text{abs}}^{\text{Cantera}}(T_{\text{ref}},Y,p_0)
$$

and then solving the Cantera HPY inversion:

$$
T = T(h_{\text{target,abs}},Y,p_0).
$$

This avoids formation-enthalpy artifacts in the current non-reacting passive-mixing model.

## Enthalpy / Species Coupling Convention

The solver uses the following convention:

$$
h \text{ is the transported thermodynamic state.}
$$

Temperature is recovered from the current transported enthalpy, current composition, and thermodynamic pressure:

$$
T = T(h,Y,p_0).
$$

Therefore, when species transport changes the composition from \(Y^n\) to \(Y^{n+1}\), the solver preserves the already-transported enthalpy field and recovers a thermodynamically consistent temperature:

$$
T^{*}
=
T(h^n,Y^{n+1},p_0).
$$

The solver must not hold the old temperature fixed and recompute enthalpy as:

$$
h^{*}
=
h(T^n,Y^{n+1},p_0).
$$

That alternative would numerically add or remove sensible enthalpy when composition changes without an energy flux or source term accounting for that change.

## Runtime Ordering

The current time step uses this ordering:

1. Update transport properties on `transport_update_interval`.
2. Advance momentum and pressure projection.
3. Advance passive species, if enabled.
4. Advance passive sensible enthalpy, if enabled.
5. Recover temperature from \(h\), current \(Y\), and \(p_0\) when Cantera thermo is enabled.
6. Refresh thermodynamic properties for output and the next step.
7. Write diagnostics and VTU/PVTU output on output steps.

In pseudocode:

```text
for step = 1 ... nsteps:

    update CFL / dt if needed

    if transport_update_interval is active:
        update mu and D_k
        keep projection density rho_f = params%rho

    advance projection step:
        compute momentum RHS
        predict u*
        compute predicted face fluxes
        solve pressure correction
        correct face fluxes and velocity

    if species enabled:
        advance passive species Y_k

    if energy enabled:
        preserve transported h
        if Cantera thermo enabled and species may have changed Y:
            sync T, cp, lambda, rho_thermo from h, current Y, p0
            restore h
        advance h with advection, conduction, and qrad/rho
        if Cantera thermo enabled:
            sync T, cp, lambda, rho_thermo from updated h and current Y
            restore h
        else:
            recover T from updated h with the constant-cp relation

    write diagnostics / output if requested
```

## Momentum Projection Method

The solver uses a fractional-step projection method.

### Momentum Predictor

A local momentum right-hand side is constructed as:

$$
\mathbf{R}_c
=
-\frac{1}{V_c}\sum_{f\in\partial c}
F_{f,c}\mathbf{u}_f^{\text{adv}}
+
\frac{1}{V_c}\sum_{f\in\partial c}
\nu_f A_f
\frac{\mathbf{u}_{\text{other}}-\mathbf{u}_c}{d_f}
-
\frac{1}{\rho_f}\nabla p_c
+
\mathbf{f}.
$$

The advected velocity \(\mathbf{u}_f^{\text{adv}}\) is selected by either the configured upwind or central scheme. The pressure gradient is evaluated by a finite-volume Gauss reconstruction:

$$
\nabla p_c
=
\frac{1}{V_c}
\sum_{f\in\partial c}
p_f \mathbf{n}_{f,c} A_f.
$$

The intermediate velocity is advanced with AB2 after the first step:

$$
\mathbf{u}_c^{*}
=
\mathbf{u}_c^n
+
\Delta t
\left(
\frac{3}{2}\mathbf{R}_c^n
-
\frac{1}{2}\mathbf{R}_c^{n-1}
\right).
$$

On the first step, where no previous right-hand side exists, the solver falls back to forward Euler:

$$
\mathbf{u}_c^{*}
=
\mathbf{u}_c^n
+
\Delta t \mathbf{R}_c^n.
$$

### Predicted Face Flux

The predicted cell-centered velocity is interpolated to faces to form a predicted volumetric flux:

$$
F_f^{*}
=
\mathbf{u}_f^{*}\cdot \mathbf{n}_f A_f.
$$

The face normal is stored in owner-to-neighbor orientation. For boundary faces, the boundary velocity condition supplies the face velocity.

For all cell-local divergence and transport operations, the face flux is reoriented so that positive flux is outward from the current cell.

### Pressure Correction Equation

The pressure correction potential is:

$$
\phi = p^{n+1} - p^n.
$$

The corrected face flux is:

$$
F_f^{n+1}
=
F_f^{*}
-
\frac{\Delta t}{\rho_f}
A_f
\frac{\phi_{\text{nb}}-\phi_{\text{owner}}}{d_f}.
$$

The pressure correction is chosen so the corrected flux divergence vanishes:

$$
\sum_{f\in\partial c} F_{f,c}^{n+1} = 0.
$$

With the solver sign convention, the discrete pressure operator is:

$$
(A\phi)_c
=
\sum_{f\in\partial c}
\frac{A_f}{d_f}
\left(\phi_c-\phi_{\text{other}}\right),
$$

and the right-hand side is:

$$
b_c
=
-\frac{\rho_f}{\Delta t}
\sum_{f\in\partial c} F_{f,c}^{*}
=
-\frac{\rho_f}{\Delta t}
V_c
(\nabla\cdot\mathbf{u}^{*})_c.
$$

The linear system is:

$$
A\phi = b.
$$

The system is solved with preconditioned conjugate gradient using a diagonal Jacobi preconditioner. Neighbor indices and geometric coefficients \(A_f/d_f\) are cached.

If a pressure Dirichlet boundary exists, it removes the null space. Otherwise, the solver pins cell 1 to remove the constant-pressure null mode.

### Velocity and Pressure Correction

After solving for \(\phi\):

$$
p_c^{n+1} = p_c^n + \phi_c,
$$

and the cell-centered velocity is corrected with a Gauss gradient of \(\phi\):

$$
\mathbf{u}_c^{n+1}
=
\mathbf{u}_c^{*}
-
\frac{\Delta t}{\rho_f}
\nabla \phi_c.
$$

The corrected conservative face fluxes are then used by species and enthalpy transport.

## Species Transport Discretization

Species are advanced explicitly using the corrected face fluxes from the projection step.

For each cell \(c\) and species \(k\), the solver accumulates:

$$
R_{Y_k,c}
=
\sum_{f\in\partial c}
\left[
-F_{f,c}Y_{k,f}^{\text{up}}
+
G_{k,f}
-
Y_{k,f}^{\text{lin}}
\sum_j G_{j,f}
\right].
$$

The first term is upwind advection:

$$
-F_{f,c}Y_{k,f}^{\text{up}}.
$$

The uncorrected diffusive contribution is:

$$
G_{k,f}
=
D_{k,f}A_f
\frac{Y_{k,\text{other}}-Y_{k,c}}{d_f}.
$$

The correction term is:

$$
-
Y_{k,f}^{\text{lin}}
\sum_j G_{j,f},
$$

where:

$$
Y_{k,f}^{\text{lin}}
=
\frac{1}{2}
\left(
Y_{k,c}+Y_{k,\text{other}}
\right).
$$

This correction enforces zero net diffusive mass flux across the species set:

$$
\sum_k
\left(
G_{k,f}
-
Y_{k,f}^{\text{lin}}
\sum_j G_{j,f}
\right)
= 0
$$

up to the consistency of the face-linear mass fractions.

The explicit update is:

$$
Y_{k,c}^{n+1}
=
Y_{k,c}^{n}
+
\Delta t
\frac{R_{Y_k,c}}{V_c}.
$$

After the update, each mass fraction is clipped:

$$
0 \le Y_{k,c}^{n+1} \le 1,
$$

and the local vector is renormalized:

$$
Y_{k,c}^{n+1}
\leftarrow
\frac{Y_{k,c}^{n+1}}
{\sum_j Y_{j,c}^{n+1}}.
$$

## Sensible-Enthalpy Transport Discretization

The transported energy variable is \(h\), not \(T\).

For each cell, the solver first stores:

$$
h_c^{\text{old}} = h_c.
$$

If Cantera thermodynamics are enabled and transported species may have changed composition, the solver synchronizes the dependent thermodynamic state from the preserved enthalpy and current composition:

$$
(T_c,c_{p,c},\lambda_c,\rho_{\text{thermo},c})
=
\operatorname{sync}(h_c^{\text{old}},Y_c,p_0).
$$

The sync is implemented as a combined Cantera recovery/property-refresh pass where possible. It then restores \(h_c = h_c^{\text{old}}\) so the thermo sync cannot alter the transported state.

For Cantera thermo runs without transported species, this pre-flux sync may be skipped after initialization because composition has not changed and the previous post-flux sync is already valid.

The finite-volume enthalpy residual is:

$$
R_{h,c}
=
\sum_{f\in\partial c}
\left[
-F_{f,c}h_f^{\text{up}}
+
\frac{\lambda_f}{\rho_f}
A_f
\frac{T_{\text{other}}-T_c}{d_f}
\right].
$$

The update is:

$$
h_c^{n+1}
=
h_c^n
+
\Delta t
\left(
\frac{R_{h,c}}{V_c}
+
\frac{q_{\text{rad},c}}{\rho_f}
\right).
$$

After the update:

1. synchronize \(T^{n+1}\), \(c_p\), \(\lambda\), and diagnostic \(\rho_{\text{thermo}}\) from \(h^{n+1}\), current \(Y\), and \(p_0\);
2. restore the transported \(h^{n+1}\) to protect it from roundoff introduced during the thermo sync.

The synchronization is logically equivalent to temperature recovery followed by property refresh, but the implementation may combine both operations into one Cantera pass and may reuse cached dependent thermo state when `h`, `Y`, and `p0` are unchanged within tight tolerances.

## Constant-\(c_p\) Fallback Thermodynamics

When Cantera thermodynamics are disabled, the solver uses:

$$
h
=
h_{\text{ref}}
+
c_p
\left(
T - T_{\text{ref}}
\right).
$$

Temperature recovery is:

$$
T
=
T_{\text{ref}}
+
\frac{h-h_{\text{ref}}}{c_p}.
$$

Thermal conductivity and heat capacity are taken from the configured constant values.

## Cantera Transport Properties

The transport-property update may use Cantera to evaluate mixture viscosity and species diffusivities. The bridge call uses:

$$
T = T_{\text{energy}}
$$

when the energy equation is enabled, otherwise the configured background temperature is used. The pressure passed to Cantera is:

$$
p_0 = p_{\text{background}}.
$$

Cantera species mass fractions are built from the solver species list. If the solver species do not sum to one, remaining mass is assigned to the bath gas, usually `N2`, when available.

The projection density remains:

$$
\rho_f = \rho_{\text{params}}.
$$

Cantera density is not fed back into the pressure projection.

## Boundary Conditions

### Velocity

Velocity boundaries support wall, symmetry, periodic, Dirichlet, and Neumann-like extrapolation behavior.

For wall or Dirichlet velocity boundaries:

$$
\mathbf{u}_b = \mathbf{u}_{\text{patch}}.
$$

For symmetry boundaries, the normal component is removed:

$$
\mathbf{u}_b
=
\mathbf{u}_{\text{int}}
-
(\mathbf{u}_{\text{int}}\cdot \mathbf{n})\mathbf{n}.
$$

Other non-Dirichlet boundaries use the interior value for face evaluation.

### Species

For species Dirichlet boundaries:

$$
Y_{k,b} = Y_{k,\text{patch}}.
$$

For non-Dirichlet species boundaries:

$$
Y_{k,b} = Y_{k,\text{int}}.
$$

### Temperature and Boundary Enthalpy

For fixed-temperature boundaries, the boundary temperature is:

$$
T_b = T_{\text{patch}}.
$$

When species and Cantera thermodynamics are enabled, the boundary enthalpy is evaluated with the boundary composition:

$$
h_b = h(T_b,Y_b,p_0).
$$

The boundary composition \(Y_b\) is obtained from the species boundary condition:

- Dirichlet species boundaries use `patch_Y`.
- Non-Dirichlet species boundaries use the adjacent interior-cell composition.

This prevents fuel or oxidizer inlets from accidentally using an interior mixture or default bath-gas mixture when converting a fixed boundary temperature to enthalpy.

For non-Cantera thermodynamics, fixed-temperature boundary enthalpy uses the constant-\(c_p\) relation:

$$
h_b
=
h_{\text{ref}}
+
c_p
\left(
T_b - T_{\text{ref}}
\right).
$$

## CFL Estimate and Time-Step Control

The cell CFL rate is computed from the absolute sum of corrected conservative face fluxes:

$$
r_c
=
\frac{1}{2V_c}
\sum_{f\in\partial c}
|F_f|.
$$

The global CFL is:

$$
\text{CFL}
=
\Delta t
\max_c r_c.
$$

When dynamic time stepping is enabled, the time step is scaled to target `max_cfl`, with a growth cap:

$$
\Delta t
\leftarrow
\min
\left(
\frac{\text{CFL}_{\max}}{\max_c r_c},
1.02\,\Delta t
\right).
$$

## Diagnostics

## Profiling Names

The terminal profiler uses inclusive timings. With nested profiling enabled, the flat rows are not additive.

Current top-level regions include:

```text
Transport_Update
Projection_Step
Species_Transport
Energy_Transport
Diagnostics_Write_Flow
Diagnostics_Write_Energy
Output_Write_VTU
```

For Cantera thermo energy runs after the combined thermo-sync optimization, the main energy children are:

```text
Energy_Cantera_PreSync
Energy_PreFlux_Exchange
Energy_Flux_Update
Energy_Cantera_PostSync
Energy_Final_Exchange
```

`Energy_Cantera_PreSync` appears only when a pre-flux thermo sync is needed, such as species-enabled runs where composition may have changed before the energy update.

The solver writes flow diagnostics including:

- maximum and RMS velocity divergence,
- net boundary flux,
- kinetic energy,
- CFL,
- pressure-iteration counts,
- maximum velocity,
- total mass based on the constant flow density,
- minimum species mass fraction.

When energy is enabled, additional diagnostics include:

- minimum, maximum, and mean temperature,
- minimum, maximum, and mean sensible enthalpy,
- minimum and maximum `qrad`,
- domain integral of `qrad`,
- maximum temperature change,
- relative enthalpy update residual.

The diagnostic mass integral uses \(\rho_f\), not \(\rho_{\text{thermo}}\).

## Current Limitations

The current staged solver is suitable for validating passive transport and the enthalpy/species coupling convention. It is not yet a full reacting low-Mach formulation.

The most important missing energy term is the species-diffusion enthalpy flux:

$$
-\nabla \cdot
\left(
\sum_k h_k \mathbf{J}_k
\right).
$$

The species solver can diffuse composition, but the energy equation does not yet transport the enthalpy carried by those diffusive species fluxes. This is acceptable for the present passive validation stage, but it must be added before treating strongly coupled multicomponent diffusion, reactions, or variable-density reacting low-Mach flow.

Do not add reactions or variable-density low-Mach coupling until the passive species and passive enthalpy tests pass with the current constant-density formulation.
