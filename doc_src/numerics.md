title: Numerical Method


# Numerical Method

This document describes the numerical method currently implemented in
LowMachReact-Hex.

It replaces older staged notes that described variable density, conservative
species/enthalpy transport, and species-enthalpy diffusion as future-only.  The
current code should be understood as:

```text
stable baseline:
  constant-density incompressible / low-Mach projection solver

active guarded path:
  non-reacting variable-density low-Mach solver using Cantera thermo density,
  conservative rho*Y and rho*h transport, and a conservative low-Mach
  divergence source

future work:
  reaction source terms, reaction heat release, physical radiation coupling,
  solver-side ideal-gas EOS, full multicomponent finite-volume diffusion, and
  fully compressible acoustic flow
```

The solver is a collocated finite-volume method on hexahedral control volumes.
The flow is advanced by a fractional-step projection method.  Species are
transported as mass fractions.  Energy is transported as sensible enthalpy.
Cantera can provide thermodynamics, transport properties, and thermodynamic
density.

Math rendering convention in this document:

```text
inline math:   \( ... \)
display math:  $$ ... $$
code symbols:  `...`
```


---

## 1. Scope of the current model

### 1.1 Implemented physics

The currently implemented and documented model includes:

```text
- transient finite-volume flow on hexahedral control volumes
- constant-density incompressible projection mode
- guarded non-reacting variable-density low-Mach projection mode
- optional passive species transport
- optional sensible-enthalpy transport
- optional species-enthalpy diffusion correction
- Fourier heat conduction driven by grad(T)
- qrad storage and source-term hook
- Cantera h <-> T thermodynamics
- Cantera cp, thermal conductivity, viscosity, diffusivity, and rho_thermo
- Cantera named-phase selection
- variable-density diagnostics for projection, continuity, and energy budgets
```

### 1.2 Not currently implemented

The current model does **not** include:

```text
- reaction source terms
- reaction heat release
- a coupled physical radiation solver
- Soret or Dufour transport
- a full coupled Stefan-Maxwell finite-volume diffusion operator
- a solver-side ideal-gas EOS path independent of Cantera
- acoustic compressibility, shocks, or fully compressible Navier-Stokes
```

`qrad` exists as a volumetric energy source hook.  It is not yet produced by a
physical radiation model unless a case or future coupling explicitly fills it.

---

## 2. Notation

For a cell \(c\):

```text
V_c       cell volume
f         face of cell c
A_f       face area
d_f       normal distance used for face gradients
n_f       stored face normal, owner-to-neighbor orientation
n_{f,c}   face normal reoriented outward from cell c
F_f       stored volumetric face flux, owner-to-neighbor orientation
F_{f,c}   volumetric face flux reoriented outward from cell c
m_dot_f   stored mass face flux, owner-to-neighbor orientation
m_dot_{f,c} mass face flux reoriented outward from cell c
```

The stored volumetric face flux is:

$$
F_f = u_f \cdot n_f A_f .
$$

The stored mass face flux is:

$$
\dot{m}_f = \rho_f F_f .
$$

During cell-local residual assembly, the solver reorients face fluxes so that
positive outward flux means flow leaving the current cell.

For an upwind transported scalar \(\psi\):

$$
\psi_f^{up} =
\begin{cases}
\psi_c, & F_{f,c} \ge 0, \\
\psi_{other}, & F_{f,c} < 0 .
\end{cases}
$$

For mass-flux-based transport in variable-density mode, the upwind direction is
based on the outward mass flux.

---

## 3. Density and pressure variables

The solver uses two different pressure concepts and two different density
concepts.

### 3.1 Projection pressure

The pressure field in the flow solver is a hydrodynamic projection pressure-like
variable.  It enforces the projection constraint:

```text
constant-density mode:
  div(u) = 0

variable-density low-Mach mode:
  div(u) = S_projection
```

This pressure is not the Cantera thermodynamic pressure.

### 3.2 Thermodynamic pressure

Cantera state calls use the uniform thermodynamic/background pressure:

$$
p_0 = \texttt{background\_press}.
$$

This pressure is used for:

```text
h(T,Y,p0)
T(h,Y,p0)
cp(T,Y,p0)
lambda(T,Y,p0)
rho(T,Y,p0)
mu(T,Y,p0)
D_k(T,Y,p0)
```

### 3.3 Flow density

The active flow density is `transport%rho`.

Constant-density mode:

$$
\rho = \rho_{params}.
$$

Variable-density mode:

$$
\texttt{transport\%rho} \leftarrow \texttt{energy\%rho\_thermo}.
$$

### 3.4 Thermodynamic density

Cantera thermo sync returns:

$$
\rho_{thermo} = \rho(T,Y,p_0).
$$

In constant-density mode, `rho_thermo` is diagnostic.  In variable-density mode,
it becomes the active flow density after the density sync.

---

## 4. Runtime modes

### 4.1 Constant-density projection mode

Input pattern:

```fortran
enable_variable_density = .false.
density_eos = "constant"
```

Constraint:

$$
\nabla \cdot u = 0 .
$$

Active density:

$$
\rho = \rho_{params}.
$$

Mass flux:

$$
\dot{m}_f = \rho_{params} F_f .
$$

Cantera may still provide thermodynamic and transport properties, but Cantera
density does not drive the projection or continuity constraint.

### 4.2 Cantera-assisted constant-density mode

Input pattern:

```fortran
enable_variable_density = .false.
density_eos = "constant"
enable_energy = .true.
enable_cantera_thermo = .true.
```

The solver still uses the constant-density projection:

$$
\nabla \cdot u = 0 .
$$

Cantera provides:

```text
T(h,Y,p0)
h(T,Y,p0)
cp
lambda
rho_thermo
mu
D_k
```

but `rho_thermo` remains diagnostic.

### 4.3 Guarded variable-density low-Mach mode

Supported input pattern:

```fortran
enable_variable_density = .true.
density_eos = "cantera"
enable_energy = .true.
enable_cantera_thermo = .true.
thermo_update_interval = 1
enable_reactions = .false.
```

Active density:

$$
\rho = \rho(T,Y,p_0)
$$

from the selected Cantera phase.

Projection target:

$$
\nabla \cdot u = S .
$$

Species transport uses a conservative \(\rho Y_k\) update.  Energy transport
uses a conservative \(\rho h\) update.  Reactions remain disabled.

The mode is called guarded because unsupported combinations are not intended to
run as production modes.

---

## 5. Momentum equation

The momentum equation is advanced in fractional-step form.

In constant-density mode, the modeled equation is:

$$
\frac{\partial u}{\partial t}
+
\nabla \cdot (u u)
=
-\frac{1}{\rho}\nabla p
+
\nabla \cdot (\nu \nabla u)
+
f .
$$

Here:

```text
rho  active flow density
nu   kinematic viscosity used by the momentum operator
p    projection pressure
f    configured body acceleration
```

In variable-density mode, the projection operator and correction use the active
density field.  The current formulation remains a low-Mach projection method,
not a fully compressible momentum equation with acoustic pressure-density
coupling.

When Cantera viscosity is active, the bridge provides dynamic viscosity \(\mu\)
and the solver may compute:

$$
\nu = \frac{\mu}{\rho}.
$$

When variable viscosity is disabled, the configured constant `nu` is preserved.

---

## 6. Momentum predictor

For cell \(c\), the momentum residual has the form:

$$
R_{u,c}
=
-\frac{1}{V_c}
\sum_{f \in \partial c}
F_{f,c} u_f^{adv}
+
\frac{1}{V_c}
\sum_{f \in \partial c}
\nu_f A_f
\frac{u_{other}-u_c}{d_f}
-
\frac{1}{\rho_c}\nabla p_c
+
f_c .
$$

The advected velocity \(u_f^{adv}\) is selected according to the configured
convection scheme, typically upwind for initial robustness or central for
cleaner laminar validation.

The pressure gradient is reconstructed by a finite-volume/Gauss form:

$$
\nabla p_c =
\frac{1}{V_c}
\sum_{f \in \partial c}
p_f n_{f,c} A_f .
$$

After the first step, the predictor uses an AB2 form:

$$
u_c^*
=
u_c^n
+
\Delta t
\left(
\frac{3}{2}R_{u,c}^n
-
\frac{1}{2}R_{u,c}^{n-1}
\right).
$$

On the first step it uses forward Euler:

$$
u_c^*
=
u_c^n + \Delta t R_{u,c}^n .
$$

---

## 7. Predicted face flux

The predicted cell velocity is interpolated to faces to form the predicted
volumetric flux:

$$
F_f^* = u_f^* \cdot n_f A_f .
$$

For boundary faces, the boundary velocity condition supplies the face velocity
where appropriate.  The stored face normal is owner-to-neighbor.  Cell-local
residuals reorient this flux outward from the current cell.

Velocity boundary conditions prescribe velocity, and therefore prescribe
volumetric flux.  In variable-density mode, the corresponding mass flux is
computed from density:

$$
\dot{m}_f = \rho_f F_f .
$$

Therefore equal inlet velocities do not generally imply equal inlet mass fluxes
when density differs.

---

## 8. Pressure projection

### 8.1 Constant-density projection

The constant-density projection corrects face fluxes as:

$$
F_f^{n+1}
=
F_f^*
-
\Delta t
\frac{A_f}{\rho d_f}
\left(
\phi_{nb} - \phi_{owner}
\right).
$$

The corrected flux must satisfy:

$$
\sum_{f \in \partial c} F_{f,c}^{n+1} = 0 .
$$

With the solver sign convention, the discrete operator is:

$$
(A\phi)_c =
\sum_{f \in \partial c}
\frac{A_f}{d_f}
\left(
\phi_c - \phi_{other}
\right),
$$

and the right-hand side is:

$$
b_c =
-\frac{\rho}{\Delta t}
\sum_{f \in \partial c} F_{f,c}^*
=
-\frac{\rho}{\Delta t}
V_c
(\nabla \cdot u^*)_c .
$$

The linear system is:

$$
A\phi = b .
$$

### 8.2 Variable-density low-Mach projection

In variable-density low-Mach mode, the target is:

$$
\nabla \cdot u = S .
$$

The corrected flux must satisfy the finite-volume constraint:

$$
\sum_{f \in \partial c} F_{f,c}^{n+1}
=
S_c V_c .
$$

The face correction uses the active face density:

$$
F_f^{n+1}
=
F_f^*
-
\Delta t
\frac{A_f}{\rho_f d_f}
\left(
\phi_{nb} - \phi_{owner}
\right).
$$

The variable-coefficient operator uses:

$$
a_f = \frac{A_f}{\rho_f d_f}.
$$

The right-hand side is:

$$
b_c =
-\frac{1}{\Delta t}
\left[
\sum_{f \in \partial c} F_{f,c}^* - S_c V_c
\right].
$$

The source used by a given projection is copied into
`fields%projection_divergence_source` immediately before projection RHS assembly.
Projection diagnostics should therefore compare the corrected fluxes against
this projection-time source, not against a source updated later in the timestep.

### 8.3 Pressure gauge

If a pressure Dirichlet boundary exists, it anchors the pressure system.

For pure-Neumann or closed/periodic constant-density cases, the constant
pressure null space is removed by the historical cell-1 pressure pin.

For the guarded variable-density path, the pure-Neumann null space is handled by
a zero-mean pressure gauge over owned cells.

### 8.4 Outlet/source compatibility

For incompressible mode with Neumann outlets, the predicted boundary volume flux
is balanced toward zero net boundary volume flux.

For variable-density low-Mach mode, the compatible target is:

$$
\int_{\partial \Omega} u \cdot n \, dA
=
\int_{\Omega} S \, dV .
$$

Outlet flux balancing therefore uses the low-Mach source integral in
variable-density mode.

---

## 9. Low-Mach divergence source

The active conservative low-Mach source in variable-density mode is:

$$
S
=
\frac{\rho^{old} - \rho}{\rho \Delta t}
-
\frac{u \cdot \nabla \rho}{\rho}.
$$

This is equivalent to targeting conservative continuity:

$$
\frac{\partial \rho}{\partial t}
+
\nabla \cdot (\rho u)
=
0 .
$$

The advective density-gradient term is evaluated explicitly from the corrected
volumetric face fluxes:

$$
(u \cdot \nabla \rho)_c
\approx
\frac{1}{V_c}
\sum_{f \in \partial c}
F_{f,c}
\left(
\rho_f - \rho_c
\right).
$$

For coupled/internal faces, \(\rho_f\) is centered/interpolated from neighboring
cell densities.  For physical boundary faces, the current staged implementation
does not yet carry a full EOS boundary-density state in this source term; the
boundary contribution uses the local cell density treatment.  A density-aware
boundary-state path is a future hardening step.

---

## 10. Species equations

### 10.1 Constant-density species equation

For each transported species \(k\):

$$
\frac{\partial Y_k}{\partial t}
+
\nabla \cdot (u Y_k)
=
\nabla \cdot (D_k \nabla Y_k),
$$

with no reaction source:

$$
\dot{\omega}_k = 0 .
$$

The finite-volume residual is assembled from corrected face fluxes:

$$
R_{Y_k,c}
=
\sum_{f \in \partial c}
\left[
-F_{f,c}Y_{k,f}^{up}
+
G_{k,f}
-
Y_{k,f}^{lin}
\sum_j G_{j,f}
\right].
$$

The uncorrected diffusive contribution is:

$$
G_{k,f}
=
D_{k,f} A_f
\frac{Y_{k,other}-Y_{k,c}}{d_f}.
$$

The correction term enforces zero net diffusive species mass flux across the
transported species set:

$$
\sum_k
\left[
G_{k,f}
-
Y_{k,f}^{lin}
\sum_j G_{j,f}
\right]
=
0
$$

up to the consistency of the face-linear mass fractions.

The explicit constant-density update is:

$$
Y_{k,c}^{n+1}
=
Y_{k,c}^{n}
+
\Delta t
\frac{R_{Y_k,c}}{V_c}.
$$

After the update, mass fractions are clipped and renormalized so the local
species vector remains bounded and sums to one.

### 10.2 Variable-density species equation

In variable-density mode, the species update uses a conservative form:

$$
\frac{\partial (\rho Y_k)}{\partial t}
+
\nabla \cdot (\rho u Y_k)
=
\nabla \cdot (\rho D_k \nabla Y_k)
$$

with the same non-reacting assumption:

$$
\dot{\omega}_k = 0 .
$$

The staged branch uses:

```text
advection: fields%mass_flux * Y_upwind
diffusion: rho_f * D_k * grad(Y_k) * area
update:    rho*Y
```

The diffusive-flux correction is applied consistently with the species
diffusion fluxes.

---

## 11. Energy equation

### 11.1 Transported variable

The transported thermodynamic variable is mixture sensible enthalpy:

$$
h \quad [\mathrm{J/kg}].
$$

Temperature is dependent:

$$
T = T(h,Y,p_0).
$$

The solver does not treat \(T\) as the transported energy state when Cantera
thermo is enabled.

### 11.2 Constant-density energy equation

In constant-density mode, the sensible enthalpy equation is interpreted as:

$$
\frac{\partial h}{\partial t}
+
\nabla \cdot (u h)
=
\frac{1}{\rho}
\nabla \cdot (\lambda \nabla T)
-
\frac{1}{\rho}
\nabla \cdot \left(\sum_k h_k J_k\right)
+
\frac{q_{rad}}{\rho}.
$$

The species-enthalpy diffusion term is included only when
`enable_species_enthalpy_diffusion = .true.`.  When that option is off, the
equation reduces to advection, conduction, and `qrad`.

Conduction is driven by temperature:

$$
q_{cond} = -\lambda \nabla T .
$$

It is not driven by \(\nabla h\).

### 11.3 Variable-density energy equation

In variable-density mode, the conservative energy update is interpreted as:

$$
\frac{\partial (\rho h)}{\partial t}
+
\nabla \cdot (\rho u h)
=
\nabla \cdot (\lambda \nabla T)
-
\nabla \cdot \left(\sum_k h_k J_k\right)
+
q_{rad}.
$$

The staged finite-volume branch uses:

```text
advection:  fields%mass_flux * h_upwind
conduction: lambda * grad(T) * area
source:     qrad
update:     rho*h
```

At the end of the update, \(h\) is recovered from the conservative state and the
active density.

---

## 12. Sensible enthalpy convention

Cantera absolute enthalpy includes formation/reference contributions.  The
solver stores sensible enthalpy relative to `energy_reference_T` at the same
composition:

$$
h_{sens}(T,Y,p_0)
=
h_{abs}^{Cantera}(T,Y,p_0)
-
h_{abs}^{Cantera}(T_{ref},Y,p_0).
$$

Temperature recovery adds back the same-composition reference enthalpy:

$$
h_{target,abs}
=
h_{sens}
+
h_{abs}^{Cantera}(T_{ref},Y,p_0),
$$

then uses Cantera's HP inversion:

$$
T = T(h_{target,abs},Y,p_0).
$$

This avoids artificial heat release during non-reacting passive mixing.

---

## 13. Enthalpy/species coupling convention

The solver uses Option A:

```text
h is the transported thermodynamic state.
T is recovered from h, Y, and p0.
```

When species transport changes composition from \(Y^n\) to \(Y^{n+1}\), the
solver preserves transported enthalpy and recovers temperature from the new
state:

$$
T^* = T(h, Y^{n+1}, p_0).
$$

The solver must not preserve old temperature by rebuilding enthalpy as:

$$
h^* = h(T^n,Y^{n+1},p_0).
$$

That alternative can add or remove sensible enthalpy when composition changes
without a corresponding energy flux or source.

---

## 14. Energy finite-volume update

Before the energy flux update, the solver stores the current transported
enthalpy.  If Cantera thermo is enabled and transported species may have changed,
it synchronizes the dependent thermodynamic state from the preserved \(h\),
current \(Y\), and \(p_0\):

$$
(T,c_p,\lambda,\rho_{thermo})
=
\operatorname{sync}(h,Y,p_0).
$$

The sync is logically equivalent to recovering \(T\) from \(h,Y,p_0\) and then
refreshing properties.  The implementation may combine both steps into one
Cantera pass.  The transported \(h\) is preserved.

The constant-density enthalpy residual has the form:

$$
R_{h,c}
=
\sum_{f \in \partial c}
\left[
-F_{f,c}h_f^{up}
+
\frac{\lambda_f}{\rho}
A_f
\frac{T_{other}-T_c}{d_f}
+
R_{hJ,f}
\right],
$$

where \(R_{hJ,f}\) represents the optional species-enthalpy diffusion
contribution when enabled.

The constant-density update is:

$$
h_c^{n+1}
=
h_c^n
+
\Delta t
\left(
\frac{R_{h,c}}{V_c}
+
\frac{q_{rad,c}}{\rho}
\right).
$$

The variable-density branch updates the conservative state \(\rho h\) using the
operator density/time level.  Direct update closure diagnostics are computed
from the same operator state.

After the update, the solver synchronizes:

$$
(T,c_p,\lambda,\rho_{thermo})
=
\operatorname{sync}(h^{n+1},Y,p_0).
$$

---

## 15. Constant-\(c_p\) fallback thermodynamics

When Cantera thermodynamics are disabled:

$$
h
=
h_{ref}
+
c_p
(T-T_{ref}).
$$

Temperature recovery is:

$$
T
=
T_{ref}
+
\frac{h-h_{ref}}{c_p}.
$$

The heat capacity and thermal conductivity are taken from the configured
constant values:

```text
energy_cp
energy_lambda
```

---

## 16. Cantera transport properties

The Cantera bridge can evaluate:

```text
mu
D_k
cp
lambda
rho_thermo
h(T,Y,p0)
T(h,Y,p0)
species sensible enthalpies
```

For transport-property updates, the temperature source is:

```text
energy enabled:
  energy%T

energy disabled:
  background_temp
```

The pressure passed to Cantera is:

$$
p_0 = \texttt{background\_press}.
$$

The bridge maps transported solver species names to Cantera species names and
constructs a full Cantera mass-fraction vector.  Missing/bath composition is
handled through the configured/default species behavior.

The finite-volume species operator currently consumes scalar \(D_k\)
coefficients.  This is a natural match for mixture-averaged diffusion.  Selecting
a Cantera multicomponent or high-pressure transport phase does not by itself
make the FV species operator a full coupled Stefan-Maxwell diffusion solver.

---

## 17. Boundary conditions

### 17.1 Velocity boundaries

Velocity boundaries support wall, moving wall, symmetry/slip, periodic,
fixed-value, and zero-gradient style behavior.

For fixed-value or wall-type velocity boundaries:

$$
u_b = u_{patch}.
$$

For symmetry boundaries, the normal component is removed:

$$
u_b =
u_{int}
-
(u_{int} \cdot n)n .
$$

Velocity boundary conditions prescribe velocity and therefore volumetric flux.
Mass flux follows from density.

### 17.2 Pressure boundaries

Pressure supports fixed-value, zero-gradient, periodic, and
symmetry-compatible behavior.

For pure-Neumann pressure systems, a gauge is needed to remove the constant
null space.

### 17.3 Species boundaries

For fixed species boundaries:

$$
Y_{k,b} = Y_{k,patch}.
$$

For non-Dirichlet species boundaries:

$$
Y_{k,b} = Y_{k,int}.
$$

### 17.4 Temperature and boundary enthalpy

For fixed-temperature boundaries:

$$
T_b = T_{patch}.
$$

When Cantera thermo and species are enabled, boundary enthalpy is evaluated with
the boundary composition:

$$
h_b = h(T_b,Y_b,p_0).
$$

Boundary composition is obtained from the species boundary condition:

```text
Dirichlet species boundary:
  Y_b = patch_Y

non-Dirichlet species boundary:
  Y_b = Y_int
```

This prevents fixed-temperature fuel/oxidizer inlets from accidentally using an
interior mixture or a default bath composition.

For constant-\(c_p\) thermodynamics:

$$
h_b =
h_{ref}
+
c_p(T_b-T_{ref}).
$$

### 17.5 Boundary density status

The current variable-density implementation uses active cell/face densities in
mass fluxes and the low-Mach source.  A fully centralized boundary EOS state,

$$
\rho_b = \rho(T_b,Y_b,p_0),
$$

is a planned hardening step for density-aware inlets, outlets, and future
mass-flow boundary conditions.

---

## 18. Radiation source convention

The volumetric energy source is:

$$
q_{rad} \quad [\mathrm{W/m^3}].
$$

Sign convention:

```text
qrad > 0: adds energy to the gas
qrad < 0: removes energy from the gas
```

If a future radiation model reports positive radiative loss, the coupling layer
should convert it as:

$$
q_{rad} = -q_{loss}.
$$

---

## 19. CFL estimate and dynamic timestep

The cell CFL rate is computed from corrected volumetric face fluxes:

$$
r_c =
\frac{1}{2V_c}
\sum_{f \in \partial c}
|F_{f,c}|.
$$

The global CFL is:

$$
\mathrm{CFL}
=
\Delta t \max_c r_c .
$$

When dynamic time stepping is enabled, the timestep targets `max_cfl` with a
growth cap:

$$
\Delta t
\leftarrow
\min\left(
\frac{\mathrm{CFL}_{max}}{\max_c r_c},
1.02\Delta t
\right).
$$

---

## 20. Main timestep ordering

The timestep is transient.  Steady solutions are obtained only when the
time-marched fields and diagnostics stop changing.

Current high-level ordering:

```text
for each timestep:

  1. Update CFL / dynamic dt if enabled.

  2. Refresh transport properties on transport_update_interval.
     - constant-density mode preserves params%rho
     - variable-density mode preserves active thermo density during transport update

  3. Advance momentum and pressure projection.
     - constant density: div(u)=0
     - variable density: div(u)=S_projection
     - corrected face_flux and mass_flux are updated

  4. Advance species, if enabled.
     - constant density: Y update
     - variable density: rho*Y update

  5. Advance sensible enthalpy, if enabled.
     - preserve transported h
     - sync thermo from h,Y,p0 when needed
     - update h or rho*h
     - sync T, cp, lambda, rho_thermo

  6. If variable-density mode is enabled:
     - sync active density from rho_thermo
     - compute conservative low-Mach source for the next projection
     - advance density history as required by the source update

  7. Write diagnostics and VTU/PVTU output on output steps.
```

Important time-level distinction:

```text
fields%projection_divergence_source:
  source copied immediately before projection RHS assembly

fields%divergence_source:
  current source after the latest density/thermo update
```

Projection validation uses the projection-time source.

---

## 21. Diagnostics

### 21.1 Constant-density projection diagnostics

In constant-density mode, raw divergence is a valid projection diagnostic:

```text
max |div(u)|
rms div(u)
net boundary volume flux
pressure iterations
pressure residual
```

### 21.2 Variable-density projection diagnostics

In variable-density low-Mach mode, raw divergence is not expected to be zero.
The target is:

$$
\nabla \cdot u = S_{projection}.
$$

Primary projection diagnostics:

```text
divu_minus_S_projection_max
divu_minus_S_projection_l2
relative_divu_minus_S_projection_max
relative_divu_minus_S_projection_l2
net_boundary_volume_flux_minus_integral_S_projection_dV
```

Current-source diagnostics remain useful for source evolution:

```text
divu_minus_S_current_*
S_current_minus_S_projection_*
net_boundary_volume_flux_minus_integral_S_current_dV
```

### 21.3 Conservative continuity diagnostics

Primary variable-density continuity diagnostics include:

```text
integral_drho_dt_plus_div_mass_flux_dV
conservative_residual_l2
relative_conservative_residual_l2
mass_balance_defect_*
```

### 21.4 Species and energy diagnostics

Species diagnostics include:

```text
sum_Y min/max/mean/l2
per-species integrals
transported_species_mass_sum
species boundary flux trends
```

Energy diagnostics include:

```text
h min/max
T min/max
rho*h integral
qrad integral
species-enthalpy diffusion integral
reported conductive boundary flux
direct energy update closure
operator-consistent rho*h budget
output-state density reconciliation
```

Primary energy closure metrics:

```text
relative_last_energy_update_balance_defect
rel_operator_recon_defect
rel_output_recon_defect
```

---

## 22. ParaView output fields

The VTU/PVTU output is for visualization and spatial debugging.  Quantitative
pass/fail decisions should use CSV diagnostics.

Common fields:

```text
velocity
pressure
rho
nu
temperature
enthalpy
qrad
cp
thermal_conductivity
rho_thermo
thermal_diffusivity
thermo_pressure
Y_<species>
D_<species>
sum_Y
mass_flux_vector
mass_flux_divergence
```

In variable-density mode, additional debug fields may include:

```text
lowmach_source_current
lowmach_source_projection
lowmach_source_difference
divu_recomputed
divu_minus_S_projection
divu_minus_S_current
rho_current
rho_projection
rho_current_minus_projection
mass_flux_divergence_recomputed
lowmach_source_history_estimate
lowmach_source_advective_density
u_dot_grad_rho
continuity_residual_estimate
rho_h_output_state
rho_h_operator_consistent
rho_h_density_reconciliation
relative_rho_h_density_reconciliation
```

Most fields are cell-centered finite-volume data.  ParaView filters such as
`Plot Over Line`, `Resample To Line`, and `Cell Data to Point Data` interpolate
or average fields.  Compare like with like when validating.

---

## 23. Current limitations and future hardening

Current limitations:

```text
- no reaction source terms
- no reaction heat release
- no physical radiation solver coupled through qrad
- no full coupled multicomponent FV diffusion operator
- no centralized EOS boundary density state for all boundary mass fluxes
- no mass-flow boundary condition yet
- scalar transport may be diffusive when using first-order/upwind paths
- low-Mach variable-density mode remains validation-driven
```

Near-term hardening targets:

```text
- density-aware boundary state evaluation: rho_b = EOS(T_b,Y_b,p0)
- optional mass-flow inlet boundary condition
- broader EOS/pressure validation matrix
- mesh/dt/MPI refinement studies for variable-density cases
- stronger scalar/enthalpy convection options after baseline validation
- restart support for long transient-to-steady runs
```

---

## 24. Summary of active equations

Constant-density flow:

$$
\nabla \cdot u = 0 .
$$

Constant-density species:

$$
\frac{\partial Y_k}{\partial t}
+
\nabla \cdot (uY_k)
=
\nabla \cdot (D_k\nabla Y_k).
$$

Constant-density sensible enthalpy:

$$
\frac{\partial h}{\partial t}
+
\nabla \cdot (u h)
=
\frac{1}{\rho}\nabla \cdot (\lambda\nabla T)
-
\frac{1}{\rho}\nabla \cdot \left(\sum_k h_k J_k\right)
+
\frac{q_{rad}}{\rho},
$$

with the species-enthalpy diffusion term present only when enabled.

Variable-density low-Mach constraint:

$$
\nabla \cdot u
=
\frac{\rho^{old} - \rho}{\rho\Delta t}
-
\frac{u \cdot \nabla \rho}{\rho}.
$$

Variable-density conservative species:

$$
\frac{\partial(\rho Y_k)}{\partial t}
+
\nabla \cdot(\rho u Y_k)
=
\nabla \cdot(\rho D_k \nabla Y_k).
$$

Variable-density conservative enthalpy:

$$
\frac{\partial(\rho h)}{\partial t}
+
\nabla \cdot(\rho u h)
=
\nabla \cdot(\lambda\nabla T)
-
\nabla \cdot \left(\sum_k h_k J_k\right)
+
q_{rad},
$$

again with the species-enthalpy diffusion term present only when enabled.

The thermodynamic state is:

$$
T = T(h,Y,p_0),
$$

with:

$$
p_0 = \texttt{background\_press}.
$$