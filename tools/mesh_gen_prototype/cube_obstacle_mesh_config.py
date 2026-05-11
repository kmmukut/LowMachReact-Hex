# Simple 3D cube-obstacle mesh.

X  = [0.0, 1.0, 1.25, 4.0]
NX = [8,   2,   22]

Y  = [0.0, 0.375, 0.625, 1.0]
NY = [3,   2,     3]

Z  = [0.0, 0.375, 0.625, 1.0]
NZ = [3,   2,     3]

# Omit central cube block.
SOLID_BLOCKS = [
    (1, 1, 1),
]

PATCHES = {
    "xmin": "inlet",
    "xmax": "outlet",
    "ymin": "wall",
    "ymax": "wall",
    "zmin": "zmin",
    "zmax": "zmax",
    "solid": "obstacle",
}

# Enable if your case.nml uses zmin/zmax periodic:
# PERIODIC = [("zmin", "zmax")]
PERIODIC = []
