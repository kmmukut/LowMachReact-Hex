#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

geo_file="${1:-counterflow.geo}"

if [[ ! -f "$geo_file" ]]; then
  echo "ERROR: geo file not found: $geo_file" >&2
  echo "Usage: $0 [geo_file]" >&2
  echo "Default: $0 counterflow.geo" >&2
  exit 1
fi

base_name="$(basename "$geo_file" .geo)"
msh_file="${base_name}.msh"

echo "Generating mesh from: $geo_file"
echo "Writing Gmsh mesh:    $msh_file"
echo "Writing native mesh:  mesh_native"

rm -rf mesh_native
gmsh -3 "$geo_file" -o "$msh_file"
../../tools/mesh/convert_gmsh_hex.py "$msh_file" mesh_native