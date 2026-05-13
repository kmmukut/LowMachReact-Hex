title: Cantera C++ Source Bridge

# Cantera C++ Source Bridge

The Cantera C++ bridge is implemented in:

- [src/cantera_interface.cpp](https://github.com/kmmukut/LowMachReact-Hex/blob/dev/src/cantera_interface.cpp)

FORD is used primarily for Fortran source documentation. The C++ bridge is linked here because FORD's non-Fortran parsing support is limited.

## Role in the solver

`cantera_interface.cpp` contains the C ABI wrappers used by the Fortran solver to call Cantera for:

- mechanism initialization
- species-name mapping
- transport-property evaluation
- sensible enthalpy evaluation
- temperature recovery from enthalpy
- combined thermo synchronization:

(T, cp, lambda, rho_thermo) = sync(h, Y, p0)

The transported energy state remains h; the C++ bridge updates only dependent thermodynamic properties.
