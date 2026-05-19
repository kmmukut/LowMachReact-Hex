#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

gmsh -3 channel.geo -o channel.msh
../../tools/mesh/convert_gmsh_hex.py channel.msh mesh_native \
  --periodic xmin:xmax \
  --periodic zmin:zmax

