# Thermo-transport demo mechanisms

These small Cantera YAML files are demo mechanisms for the named-phase selector
added by Patch 072.

Each file contains one phase named:

```text
gas
```

The filename identifies the thermo/EOS and transport model:

```text
ideal_mixavg.yaml   ideal-gas       + mixture-averaged transport
ideal_multi.yaml    ideal-gas       + multicomponent transport
pr_high.yaml        Peng-Robinson   + high-pressure transport
pr_chung.yaml       Peng-Robinson   + high-pressure-Chung transport
rk_high.yaml        Redlich-Kwong   + high-pressure transport
rk_chung.yaml       Redlich-Kwong   + high-pressure-Chung transport
```

All demo files use the same small stable species set imported from `gri30.yaml`:

```text
CH4, O2, N2, CO2
```

The YAML `state` block is only Cantera's default load state. The solver later
sets cell states with its own temperature, pressure, and mass fractions.

Use from `case.nml` like:

```fortran
&fluid_input
  density_eos = "cantera"
  cantera_mech_file = "mechanisms/pr_high.yaml"
  cantera_phase_name = "gas"
/
```

Before using a demo mechanism in the solver, test it directly:

```bash
python - <<'PY'
import cantera as ct
for fname in [
    "mechanisms/ideal_mixavg.yaml",
    "mechanisms/ideal_multi.yaml",
    "mechanisms/pr_high.yaml",
    "mechanisms/pr_chung.yaml",
    "mechanisms/rk_high.yaml",
    "mechanisms/rk_chung.yaml",
]:
    print("\n---", fname)
    gas = ct.Solution(fname, "gas")
    gas.TPX = 300.0, ct.one_atm, {"CH4": 0.05, "O2": 0.21, "N2": 0.74}
    print("thermo:", gas.thermo_model)
    print("transport:", gas.transport_model)
    print("rho:", gas.density)
    print("cp:", gas.cp_mass)
    print("lambda:", gas.thermal_conductivity)
    print("mu:", gas.viscosity)
PY
```

For detailed mechanisms created from CHEMKIN, use Cantera's converter, for
example:

```bash
ck2yaml --input chem.inp --thermo therm.dat --transport tran.dat --output mechanism.yaml
```

See `doc_src/input_configuration_guide.md` for more details.
