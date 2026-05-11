# Optimized 2D square-obstacle mesh.
# One-cell-thick in z. Good for your 2D-style flow-over-square case.

X  = [0.0, 0.75, 1.0, 1.25, 2.0, 4.0]
NX = [12,  8,    8,   24,   32]

Y  = [0.0, 0.25, 0.375, 0.625, 0.75, 1.0]
NY = [4,   4,    8,     4,     4]

Z  = [0.0, 1.0 / 32.0]
NZ = [1]

# Omit the block x=1.0..1.25, y=0.375..0.625, z=0..1/32.
# Block indices are zero-based: [i, j, k].
SOLID_BLOCKS = [
    (2, 2, 0),
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

# For this 2D-style case, use symmetric zmin/zmax in case.nml, not periodic.
PERIODIC = []
