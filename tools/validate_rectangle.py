import glob
import re
import numpy as np
import pandas as pd
import pyvista as pv
from scipy.signal import find_peaks
import matplotlib.pyplot as plt


# ----------------------------
# User settings
# ----------------------------
vtu_pattern = "./cases/rectangle_2D/output/*.vtu"

rho = 1.0
nu = 1.0e-2
mu = rho * nu

Uref = 10.0
D = 0.25
Lz = 1.0 / 32.0
Aref = D * Lz

dt_solver = 1.0e-4
output_interval = 100
dt_output = dt_solver * output_interval

# obstacle limits
x1, x2 = 1.0, 1.25
y1, y2 = 0.375, 0.625

# tolerance for identifying obstacle surface
tol = 1.0e-8


# ----------------------------
# Helper functions
# ----------------------------
def natural_sort_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def get_array_name(mesh, candidates):
    names = list(mesh.point_data.keys()) + list(mesh.cell_data.keys())
    for c in candidates:
        if c in names:
            return c
    raise KeyError(f"Could not find any of {candidates}. Available arrays: {names}")


def cell_centers_array(surface):
    centers = surface.cell_centers().points
    return centers[:, 0], centers[:, 1], centers[:, 2]


def obstacle_face_mask(surface):
    cx, cy, cz = cell_centers_array(surface)

    on_left = (
        np.isclose(cx, x1, atol=tol)
        & (cy >= y1 - tol)
        & (cy <= y2 + tol)
    )

    on_right = (
        np.isclose(cx, x2, atol=tol)
        & (cy >= y1 - tol)
        & (cy <= y2 + tol)
    )

    on_bottom = (
        np.isclose(cy, y1, atol=tol)
        & (cx >= x1 - tol)
        & (cx <= x2 + tol)
    )

    on_top = (
        np.isclose(cy, y2, atol=tol)
        & (cx >= x1 - tol)
        & (cx <= x2 + tol)
    )

    return on_left | on_right | on_bottom | on_top


def compute_force_from_vtu(filename):
    mesh = pv.read(filename)

    # Try to detect common field names
    p_name = get_array_name(mesh, ["p", "pressure", "Pressure"])
    u_name = get_array_name(mesh, ["U", "velocity", "Velocity", "u"])

    # Surface extraction
    surf = mesh.extract_surface()
    surf = surf.compute_normals(
        cell_normals=True,
        point_normals=False,
        auto_orient_normals=True,
        consistent_normals=True,
    )

    mask = obstacle_face_mask(surf)
    obs = surf.extract_cells(mask)

    if obs.n_cells == 0:
        raise RuntimeError(f"No obstacle faces found in {filename}")

    # Cell areas
    obs = obs.compute_cell_sizes(length=False, area=True, volume=False)
    area = obs.cell_data["Area"]

    # Cell normals
    normals = obs.cell_data["Normals"]

    # Pressure field
    if p_name in obs.point_data:
        p_cell = obs.point_data_to_cell_data()[p_name]
    else:
        p_cell = obs.cell_data[p_name]

    # Pressure force: F = integral(-p n dA)
    Fp = np.sum((-p_cell[:, None] * normals) * area[:, None], axis=0)

    # If you do not have velocity gradients, this script computes pressure force only.
    # Viscous force requires grad(U) near the wall.
    F_total = Fp

    Fx = F_total[0]
    Fy = F_total[1]

    q = 0.5 * rho * Uref**2

    Cd = Fx / (q * Aref)
    Cl = Fy / (q * Aref)

    return Fx, Fy, Cd, Cl


# ----------------------------
# Process all VTU files
# ----------------------------
files = sorted(glob.glob(vtu_pattern), key=natural_sort_key)

rows = []

for i, f in enumerate(files):
    time = i * dt_output
    Fx, Fy, Cd, Cl = compute_force_from_vtu(f)

    rows.append(
        {
            "file": f,
            "time": time,
            "Fx": Fx,
            "Fy": Fy,
            "Cd": Cd,
            "Cl": Cl,
        }
    )

df = pd.DataFrame(rows)
df.to_csv("force_coefficients.csv", index=False)

print(df.head())
print(df.tail())