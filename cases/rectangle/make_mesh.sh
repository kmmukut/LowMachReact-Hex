#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

gmsh -3 rectangle.geo -format msh2 -o rectangle.msh

../../tools/mesh/convert_gmsh_hex.py \
  rectangle.msh \
  mesh_native 
  # --periodic zmin:zmax