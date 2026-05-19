#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

gmsh -3 cavity.geo -o cavity.msh
../../tools/mesh/convert_gmsh_hex.py cavity.msh mesh_native

