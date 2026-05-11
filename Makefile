# ============================================================================
# LowMachReact-Hex Makefile
#
# Important runtime issue:
# ------------------------
# Cantera from Conda depends on Conda's newer libstdc++.so.6.
# On this cluster, LD_LIBRARY_PATH can cause the loader to pick the cluster GCC
# libstdc++ first, which breaks with:
#
#   GLIBCXX_3.4.32 not found
#   CXXABI_1.3.15 not found
#
# Therefore:
#   - The real compiled ELF binary is build/lowmach_react_hex.bin
#   - The top-level ./lowmach_react_hex is a small launcher script
#   - The launcher prepends $CONDA_PREFIX/lib to LD_LIBRARY_PATH automatically
#
# Users still run:
#
#   ./lowmach_react_hex cases/channel_flow/case.nml
#
# or:
#
#   mpirun -np 4 ./lowmach_react_hex cases/channel_flow/case.nml
#
# No manual LD_LIBRARY_PATH command is needed.
# ============================================================================

SHELL := /bin/bash

EXEC      := lowmach_react_hex
BUILD_DIR := build
SRC_DIR   := src
BIN       := $(BUILD_DIR)/$(EXEC).bin

MPIFC ?= mpifort
BUILD ?= debug
NP    ?= 1
CASES_DIR ?= cases

# Automatically detect cases that have cases/<case_name>/case.nml
CASE_NAMES := $(patsubst $(CASES_DIR)/%/case.nml,%,$(wildcard $(CASES_DIR)/*/case.nml))
CASE_RELEASE_NAMES := $(addsuffix -release,$(CASE_NAMES))

# ----------------------------------------------------------------------------
# Conda / Cantera settings
# ----------------------------------------------------------------------------
#
# CONDA_PREFIX must come from:
#
#   conda activate react_env
#
# The launcher script stores this build-time path as a fallback, but at runtime
# it prefers the active CONDA_PREFIX if one is set.

CONDA_LIB := $(CONDA_PREFIX)/lib
BUILD_CONDA_PREFIX := $(CONDA_PREFIX)

# Prefer the Conda C++ compiler if it exists. Otherwise fall back to system g++.
# This is only used for cantera_interface.cpp.
CONDA_CXX := $(firstword $(wildcard \
  $(CONDA_PREFIX)/bin/x86_64-conda-linux-gnu-c++ \
  $(CONDA_PREFIX)/bin/x86_64-conda-linux-gnu-g++))

ifeq ($(CONDA_CXX),)
  CXX ?= g++
else
  CXX ?= $(CONDA_CXX)
endif

# Get Cantera include and link flags from pkg-config.
# conda activate react_env should make this pkg-config entry visible.
CANTERA_CXXFLAGS := $(shell pkg-config --cflags cantera 2>/dev/null)
CANTERA_LIBS     := $(shell pkg-config --libs cantera 2>/dev/null)

# Try to collect MPI library directories for RPATH. This works for OpenMPI.
# If the MPI wrapper does not support --showme:link, this becomes empty.
MPI_LIB_DIRS := $(shell $(MPIFC) --showme:link 2>/dev/null | tr ' ' '\n' | sed -n 's/^-L//p' | paste -sd: -)

ifneq ($(MPI_LIB_DIRS),)
  FINAL_RPATH := $(CONDA_LIB):$(MPI_LIB_DIRS)
else
  FINAL_RPATH := $(CONDA_LIB)
endif

# ----------------------------------------------------------------------------
# Compiler flags
# ----------------------------------------------------------------------------

FFLAGS_COMMON  := -std=f2008 -Wall -Wextra -Wtabs -ffree-line-length-none -J$(BUILD_DIR)
FFLAGS_DEBUG   := -O0 -g -fcheck=all -fbacktrace -pedantic
FFLAGS_RELEASE := -O3 -march=native -fno-omit-frame-pointer

CXXFLAGS := -O2 -g -Wall $(CANTERA_CXXFLAGS)

ifeq ($(BUILD),release)
  FFLAGS := $(FFLAGS_COMMON) $(FFLAGS_RELEASE)
else
  FFLAGS := $(FFLAGS_COMMON) $(FFLAGS_DEBUG)
endif

# ----------------------------------------------------------------------------
# Link flags
# ----------------------------------------------------------------------------
#
# Notes:
#
# 1. -L$(CONDA_LIB) is placed first so the linker sees Conda libraries first.
#
# 2. -Wl,--no-as-needed -lstdc++ forces libstdc++.so.6 to appear as a direct
#    NEEDED dependency of the executable. This helps ensure the correct C++
#    runtime is loaded early.
#
# 3. We still add RPATH and then patch it after linking. The launcher is the
#    most reliable fix, but RPATH is useful for diagnostics and for environments
#    without hostile LD_LIBRARY_PATH settings.

LDFLAGS := \
  -L$(CONDA_LIB) \
  -Wl,--disable-new-dtags \
  -Wl,-rpath,$(CONDA_LIB) \
  -Wl,-rpath-link,$(CONDA_LIB) \
  $(CANTERA_LIBS) \
  -Wl,--no-as-needed -lstdc++ -Wl,--as-needed \
  -lm -ldl -lpthread

# ----------------------------------------------------------------------------
# Sources and objects
# ----------------------------------------------------------------------------

F_SRCS := \
  $(SRC_DIR)/mod_kinds.f90 \
  $(SRC_DIR)/mod_input.f90 \
  $(SRC_DIR)/mod_mesh_types.f90 \
  $(SRC_DIR)/mod_mesh_io.f90 \
  $(SRC_DIR)/mod_mpi_flow.f90 \
  $(SRC_DIR)/mod_mpi_radiation.f90 \
  $(SRC_DIR)/mod_bc.f90 \
  $(SRC_DIR)/mod_fields.f90 \
  $(SRC_DIR)/mod_flow_projection.f90 \
  $(SRC_DIR)/mod_transport_properties.f90 \
  $(SRC_DIR)/mod_species.f90 \
  $(SRC_DIR)/mod_output.f90 \
  $(SRC_DIR)/main.f90

CXX_SRCS := \
  $(SRC_DIR)/cantera_interface.cpp

F_OBJS   := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(F_SRCS))
CXX_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(CXX_SRCS))
OBJS     := $(F_OBJS) $(CXX_OBJS)

.PHONY: all clean debug release check-build-env check-runtime FORCE \
        mesh-cavity mesh-channel cavity channel cavity-release channel-release \
        list-cases \
        $(CASE_NAMES) $(CASE_RELEASE_NAMES)

# ============================================================================
# Build rules
# ============================================================================

all: check-build-env $(EXEC)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Build the real executable as build/lowmach_react_hex.bin.
# Then use patchelf to force an old-style RPATH.
#
# The launcher script is still the primary fix for your cluster because
# LD_LIBRARY_PATH from modules can otherwise override RUNPATH behavior.
$(BIN): $(BUILD_DIR) $(OBJS)
	$(MPIFC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS)
	patchelf --force-rpath --set-rpath "$(FINAL_RPATH)" $@
	@echo "Dynamic path tags in $@:"
	@readelf -d $@ | grep -E 'RPATH|RUNPATH' || true

# Generate the user-facing launcher.
#
# Users run ./lowmach_react_hex exactly as before.
# The launcher automatically puts Conda's lib directory first, then execs the
# real ELF binary.
$(EXEC): $(BIN)
	@printf '%s\n' \
	'#!/usr/bin/env bash' \
	'set -euo pipefail' \
	'' \
	'# Auto-generated by Makefile. Do not edit by hand.' \
	'# This wrapper makes Conda runtime libraries take priority over cluster GCC libraries.' \
	'' \
	'BUILD_CONDA_PREFIX="$(BUILD_CONDA_PREFIX)"' \
	'RUNTIME_CONDA_PREFIX="$${CONDA_PREFIX:-$$BUILD_CONDA_PREFIX}"' \
	'' \
	'if [[ -z "$$RUNTIME_CONDA_PREFIX" || ! -d "$$RUNTIME_CONDA_PREFIX/lib" ]]; then' \
	'  echo "ERROR: Conda library directory not found: $$RUNTIME_CONDA_PREFIX/lib" >&2' \
	'  echo "Activate the Conda environment used for Cantera, then run again." >&2' \
	'  exit 1' \
	'fi' \
	'' \
	'SCRIPT_DIR="$$(cd "$$(dirname "$${BASH_SOURCE[0]}")" && pwd)"' \
	'export LD_LIBRARY_PATH="$$RUNTIME_CONDA_PREFIX/lib:$${LD_LIBRARY_PATH:-}"' \
	'' \
	'exec "$$SCRIPT_DIR/$(BIN)" "$$@"' \
	> $@
	@chmod +x $@
	@echo "Created launcher: ./$@"

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(MPIFC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ============================================================================
# Fortran module dependencies
# ============================================================================

$(BUILD_DIR)/mod_input.o: $(BUILD_DIR)/mod_kinds.o
$(BUILD_DIR)/mod_mesh_types.o: $(BUILD_DIR)/mod_kinds.o
$(BUILD_DIR)/mod_mesh_io.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o
$(BUILD_DIR)/mod_mpi_flow.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o
$(BUILD_DIR)/mod_mpi_radiation.o: $(BUILD_DIR)/mod_kinds.o
$(BUILD_DIR)/mod_bc.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mesh_types.o
$(BUILD_DIR)/mod_fields.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_bc.o
$(BUILD_DIR)/mod_flow_projection.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_bc.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_input.o
$(BUILD_DIR)/mod_transport_properties.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/cantera_interface.o
$(BUILD_DIR)/mod_species.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_bc.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_transport_properties.o
$(BUILD_DIR)/mod_output.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_species.o
$(BUILD_DIR)/main.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mesh_io.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_mpi_radiation.o $(BUILD_DIR)/mod_bc.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_transport_properties.o $(BUILD_DIR)/mod_species.o $(BUILD_DIR)/mod_output.o

# ============================================================================
# Convenience build targets
# ============================================================================

debug:
	$(MAKE) clean
	$(MAKE) BUILD=debug all

release:
	$(MAKE) clean
	$(MAKE) BUILD=release all

clean:
	rm -rf $(BUILD_DIR) $(EXEC)

FORCE:

# ============================================================================
# Environment checks
# ============================================================================

check-build-env:
	@test -n "$(CONDA_PREFIX)" || \
	  (echo "ERROR: CONDA_PREFIX is not set. Run: conda activate react_env"; exit 1)
	@test -d "$(CONDA_LIB)" || \
	  (echo "ERROR: Conda lib directory not found: $(CONDA_LIB)"; exit 1)
	@command -v $(MPIFC) >/dev/null 2>&1 || \
	  (echo "ERROR: $(MPIFC) not found."; exit 1)
	@command -v $(CXX) >/dev/null 2>&1 || \
	  (echo "ERROR: C++ compiler not found: $(CXX)"; exit 1)
	@command -v pkg-config >/dev/null 2>&1 || \
	  (echo "ERROR: pkg-config not found. Install it in react_env."; exit 1)
	@pkg-config --exists cantera || \
	  (echo "ERROR: pkg-config cannot find Cantera. Is react_env activated?"; exit 1)
	@command -v patchelf >/dev/null 2>&1 || \
	  (echo "ERROR: patchelf not found. Install with: conda install -c conda-forge patchelf"; exit 1)
	@test -f "$(CONDA_LIB)/libstdc++.so.6" || \
	  (echo "ERROR: $(CONDA_LIB)/libstdc++.so.6 not found."; exit 1)

# Check what the wrapper will make the dynamic loader choose.
check-runtime: $(EXEC)
	@echo "RPATH/RUNPATH for real binary:"
	@readelf -d $(BIN) | grep -E 'RPATH|RUNPATH' || true
	@echo
	@echo "libstdc++ selected when using launcher environment:"
	@env LD_LIBRARY_PATH="$(CONDA_LIB):$${LD_LIBRARY_PATH:-}" ldd $(BIN) | grep 'libstdc++' || true
	@echo
	@echo "Cantera selected when using launcher environment:"
	@env LD_LIBRARY_PATH="$(CONDA_LIB):$${LD_LIBRARY_PATH:-}" ldd $(BIN) | grep 'libcantera' || true

# ============================================================================
# Mesh generation
# ============================================================================

check-mesh-env:
	@command -v gmsh >/dev/null 2>&1 || \
	  (echo "ERROR: gmsh not found. Activate react_env."; exit 1)
	@python -c "import meshio" >/dev/null 2>&1 || \
	  (echo "ERROR: meshio not found. Activate react_env."; exit 1)

mesh-cavity: check-mesh-env
	./cases/lid_driven_cavity/make_mesh.sh

mesh-channel: check-mesh-env
	./cases/channel_flow/make_mesh.sh

# ============================================================================
# Convenience aliases for old case names
# ============================================================================

cavity:
	$(MAKE) lid_driven_cavity NP=$(NP)

channel:
	$(MAKE) channel_flow NP=$(NP)

cavity-release:
	$(MAKE) lid_driven_cavity-release NP=$(NP)

channel-release:
	$(MAKE) channel_flow-release NP=$(NP)

# ============================================================================
# Generic case runner
# ============================================================================

$(CASE_NAMES): FORCE
	$(MAKE) debug
	@echo "Running debug case: $(CASES_DIR)/$@/case.nml"
	mpirun -np $(NP) ./$(EXEC) "$(CASES_DIR)/$@/case.nml"

$(CASE_RELEASE_NAMES): %-release: FORCE
	$(MAKE) release
	@echo "Running release case: $(CASES_DIR)/$*/case.nml"
	mpirun -np $(NP) ./$(EXEC) "$(CASES_DIR)/$*/case.nml"

list-cases:
	@echo "Available cases:"
	@for c in $(CASE_NAMES); do echo "  $$c"; done