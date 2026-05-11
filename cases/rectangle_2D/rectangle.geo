// ============================================================================
// Optimized 2D flow-over-square-obstacle mesh
//
// One-cell-thick 3D mesh for 2D-style FV solver.
// Fine near square and wake, coarser far away.
// Axis-aligned conformal hex mesh.
// ============================================================================

SetFactory("Built-in");

Mesh.MshFileVersion = 2.2;

// ---------------------------------------------------------------------------
// Domain
// ---------------------------------------------------------------------------

Lx = 4.0;
Ly = 1.0;

// Fine near-obstacle spacing
hf = 1.0 / 32.0;

// One-cell-thick 2D-style extrusion
Lz = hf;
nz = 1;

lc = hf;

// ---------------------------------------------------------------------------
// Square obstacle
// ---------------------------------------------------------------------------

x1 = 1.0;
x2 = 1.25;

y1 = 0.375;
y2 = 0.625;

// ---------------------------------------------------------------------------
// Block coordinates
//
// x blocks:
//   [0,    0.75]  coarse upstream
//   [0.75, 1.0 ]  fine upstream of obstacle
//   [1.0,  1.25]  obstacle-width block, omitted in center row
//   [1.25, 2.0 ]  fine near wake
//   [2.0,  4.0 ]  coarse far wake
//
// y blocks:
//   [0,     0.25 ] coarse lower far field
//   [0.25,  0.375] fine lower obstacle approach
//   [0.375, 0.625] obstacle-height block, omitted in center column
//   [0.625, 0.75 ] fine upper obstacle approach
//   [0.75,  1.0  ] coarse upper far field
// ---------------------------------------------------------------------------

X[] = {0.0, 0.75, x1, x2, 2.0, Lx};
Y[] = {0.0, 0.25, y1, y2, 0.75, Ly};

// Cell counts per block.
// Fine blocks have dx = dy = 1/32.
// Coarse far blocks have dx or dy = 1/16.
NXC[] = {12, 8, 8, 24, 32};
NYC[] = {4, 4, 8, 4, 4};

NxPts = 6;
NyPts = 6;
NxSeg = 5;
NySeg = 5;

// ---------------------------------------------------------------------------
// Points
//
// point id = 1 + i + NxPts*j
// ---------------------------------------------------------------------------

For j In {0:5}
  For i In {0:5}
    p = 1 + i + NxPts*j;
    Point(p) = {X[i], Y[j], 0.0, lc};
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Lines in x direction
//
// x-line id = 1000 + i + NxSeg*j
// ---------------------------------------------------------------------------

For j In {0:5}
  For i In {0:4}
    l = 1000 + i + NxSeg*j;
    p0 = 1 + i     + NxPts*j;
    p1 = 1 + i + 1 + NxPts*j;
    Line(l) = {p0, p1};
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Lines in y direction
//
// y-line id = 2000 + i + NxPts*j
// ---------------------------------------------------------------------------

For j In {0:4}
  For i In {0:5}
    l = 2000 + i + NxPts*j;
    p0 = 1 + i + NxPts*j;
    p1 = 1 + i + NxPts*(j + 1);
    Line(l) = {p0, p1};
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Transfinite curves
// ---------------------------------------------------------------------------

// Horizontal curves
For j In {0:5}
  For i In {0:4}
    l = 1000 + i + NxSeg*j;
    Transfinite Curve {l} = NXC[i] + 1;
  EndFor
EndFor

// Vertical curves
For j In {0:4}
  For i In {0:5}
    l = 2000 + i + NxPts*j;
    Transfinite Curve {l} = NYC[j] + 1;
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Fluid surfaces
//
// Logical block i=2, j=2 is omitted.
// That omitted block is the solid square obstacle.
// ---------------------------------------------------------------------------

fluid_surfaces[] = {};

For j In {0:4}
  For i In {0:4}

    If (i != 2 || j != 2)

      s = 3000 + i + NxSeg*j;

      p00 = 1 + i     + NxPts*j;
      p10 = 1 + i + 1 + NxPts*j;
      p11 = 1 + i + 1 + NxPts*(j + 1);
      p01 = 1 + i     + NxPts*(j + 1);

      lx0 = 1000 + i + NxSeg*j;
      ly1 = 2000 + i + 1 + NxPts*j;
      lx1 = 1000 + i + NxSeg*(j + 1);
      ly0 = 2000 + i + NxPts*j;

      Curve Loop(s) = {lx0, ly1, -lx1, -ly0};
      Plane Surface(s) = {s};

      Transfinite Surface {s} = {p00, p10, p11, p01};
      Recombine Surface {s};

      fluid_surfaces[] += {s};

    EndIf

  EndFor
EndFor

// ---------------------------------------------------------------------------
// Extrude to one-cell-thick 3D mesh
// ---------------------------------------------------------------------------

Extrude {0, 0, Lz} {
  Surface{fluid_surfaces[]};
  Layers{nz};
  Recombine;
}

// ---------------------------------------------------------------------------
// Physical patches
// ---------------------------------------------------------------------------

eps = 1.0e-8;

inlet() = Surface In BoundingBox {
  -eps, -eps, -eps,
   eps, Ly + eps, Lz + eps
};

outlet() = Surface In BoundingBox {
  Lx - eps, -eps, -eps,
  Lx + eps, Ly + eps, Lz + eps
};

zmin() = Surface In BoundingBox {
  -eps, -eps, -eps,
  Lx + eps, Ly + eps, eps
};

zmax() = Surface In BoundingBox {
  -eps, -eps, Lz - eps,
  Lx + eps, Ly + eps, Lz + eps
};

bottom_wall() = Surface In BoundingBox {
  -eps, -eps, -eps,
  Lx + eps, eps, Lz + eps
};

top_wall() = Surface In BoundingBox {
  -eps, Ly - eps, -eps,
  Lx + eps, Ly + eps, Lz + eps
};

obs_bottom() = Surface In BoundingBox {
  x1 - eps, y1 - eps, -eps,
  x2 + eps, y1 + eps, Lz + eps
};

obs_top() = Surface In BoundingBox {
  x1 - eps, y2 - eps, -eps,
  x2 + eps, y2 + eps, Lz + eps
};

obs_left() = Surface In BoundingBox {
  x1 - eps, y1 - eps, -eps,
  x1 + eps, y2 + eps, Lz + eps
};

obs_right() = Surface In BoundingBox {
  x2 - eps, y1 - eps, -eps,
  x2 + eps, y2 + eps, Lz + eps
};

fluid() = Volume In BoundingBox {
  -eps, -eps, -eps,
  Lx + eps, Ly + eps, Lz + eps
};

Physical Surface("inlet") = {inlet()};
Physical Surface("outlet") = {outlet()};

Physical Surface("wall") = {
  bottom_wall(),
  top_wall(),
  obs_bottom(),
  obs_top(),
  obs_left(),
  obs_right()
};

Physical Surface("zmin") = {zmin()};
Physical Surface("zmax") = {zmax()};

Physical Volume("fluid") = {fluid()};