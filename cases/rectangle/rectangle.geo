// ============================================================================
// 3D flow-over-cube-obstacle test with perfect cube cells
// ============================================================================

SetFactory("Built-in");

Mesh.MshFileVersion = 2.2;

h = 1.0 / 8.0;

// Domain
Lx = 4.0;
Ly = 1.0;
Lz = 1.0;

// Cube obstacle aligned to grid
x1 = 1.0;     // 8h
x2 = 1.25;    // 10h, cube width = 2h

y1 = 0.375;   // 3h
y2 = 0.625;   // 5h, cube height = 2h

z1 = 0.375;   // 3h
z2 = 0.625;   // 5h, cube depth = 2h

// Mesh counts
nx1 = 8;
nx2 = 2;
nx3 = 22;

ny1 = 3;
ny2 = 2;
ny3 = 3;

nz1 = 3;
nz2 = 2;
nz3 = 3;

lc = h;

// Coordinate planes
X[] = {0, x1, x2, Lx};
Y[] = {0, y1, y2, Ly};
Z[] = {0, z1, z2, Lz};

// ---------------------------------------------------------------------------
// Points
// Point id = 1 + i + 4*j + 16*k
// ---------------------------------------------------------------------------

For k In {0:3}
  For j In {0:3}
    For i In {0:3}
      p = 1 + i + 4*j + 16*k;
      Point(p) = {X[i], Y[j], Z[k], lc};
    EndFor
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Lines in x direction
// lx id = 1000 + i + 3*j + 12*k
// ---------------------------------------------------------------------------

For k In {0:3}
  For j In {0:3}
    For i In {0:2}
      l = 1000 + i + 3*j + 12*k;
      p0 = 1 + i     + 4*j + 16*k;
      p1 = 1 + i + 1 + 4*j + 16*k;
      Line(l) = {p0, p1};
    EndFor
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Lines in y direction
// ly id = 2000 + i + 4*j + 12*k
// ---------------------------------------------------------------------------

For k In {0:3}
  For j In {0:2}
    For i In {0:3}
      l = 2000 + i + 4*j + 12*k;
      p0 = 1 + i + 4*j     + 16*k;
      p1 = 1 + i + 4*(j+1) + 16*k;
      Line(l) = {p0, p1};
    EndFor
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Lines in z direction
// lz id = 3000 + i + 4*j + 16*k
// ---------------------------------------------------------------------------

For k In {0:2}
  For j In {0:3}
    For i In {0:3}
      l = 3000 + i + 4*j + 16*k;
      p0 = 1 + i + 4*j + 16*k;
      p1 = 1 + i + 4*j + 16*(k+1);
      Line(l) = {p0, p1};
    EndFor
  EndFor
EndFor

// ---------------------------------------------------------------------------
// xy surfaces, normal z
// sxy id = 4000 + i + 3*j + 9*k
// ---------------------------------------------------------------------------

For k In {0:3}
  For j In {0:2}
    For i In {0:2}
      s = 4000 + i + 3*j + 9*k;

      lx0 = 1000 + i + 3*j     + 12*k;
      ly1 = 2000 + i + 1 + 4*j + 12*k;
      lx1 = 1000 + i + 3*(j+1) + 12*k;
      ly0 = 2000 + i + 4*j     + 12*k;

      Curve Loop(s) = {lx0, ly1, -lx1, -ly0};
      Plane Surface(s) = {s};
    EndFor
  EndFor
EndFor

// ---------------------------------------------------------------------------
// xz surfaces, normal y
// sxz id = 5000 + i + 3*k + 9*j
// ---------------------------------------------------------------------------

For j In {0:3}
  For k In {0:2}
    For i In {0:2}
      s = 5000 + i + 3*k + 9*j;

      lx0 = 1000 + i     + 3*j + 12*k;
      lz1 = 3000 + i + 1 + 4*j + 16*k;
      lx1 = 1000 + i     + 3*j + 12*(k+1);
      lz0 = 3000 + i     + 4*j + 16*k;

      Curve Loop(s) = {lx0, lz1, -lx1, -lz0};
      Plane Surface(s) = {s};
    EndFor
  EndFor
EndFor

// ---------------------------------------------------------------------------
// yz surfaces, normal x
// syz id = 6000 + j + 3*k + 9*i
// ---------------------------------------------------------------------------

For i In {0:3}
  For k In {0:2}
    For j In {0:2}
      s = 6000 + j + 3*k + 9*i;

      ly0 = 2000 + i + 4*j     + 12*k;
      lz1 = 3000 + i + 4*(j+1) + 16*k;
      ly1 = 2000 + i + 4*j     + 12*(k+1);
      lz0 = 3000 + i + 4*j     + 16*k;

      Curve Loop(s) = {ly0, lz1, -ly1, -lz0};
      Plane Surface(s) = {s};
    EndFor
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Transfinite curves
// ---------------------------------------------------------------------------

// x-direction curves
For k In {0:3}
  For j In {0:3}
    Transfinite Curve {1000 + 0 + 3*j + 12*k} = nx1 + 1;
    Transfinite Curve {1000 + 1 + 3*j + 12*k} = nx2 + 1;
    Transfinite Curve {1000 + 2 + 3*j + 12*k} = nx3 + 1;
  EndFor
EndFor

// y-direction curves
For k In {0:3}
  For i In {0:3}
    Transfinite Curve {2000 + i + 4*0 + 12*k} = ny1 + 1;
    Transfinite Curve {2000 + i + 4*1 + 12*k} = ny2 + 1;
    Transfinite Curve {2000 + i + 4*2 + 12*k} = ny3 + 1;
  EndFor
EndFor

// z-direction curves
For j In {0:3}
  For i In {0:3}
    Transfinite Curve {3000 + i + 4*j + 16*0} = nz1 + 1;
    Transfinite Curve {3000 + i + 4*j + 16*1} = nz2 + 1;
    Transfinite Curve {3000 + i + 4*j + 16*2} = nz3 + 1;
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Transfinite and recombine all surfaces
// ---------------------------------------------------------------------------

For k In {0:3}
  For j In {0:2}
    For i In {0:2}
      s = 4000 + i + 3*j + 9*k;
      Transfinite Surface {s};
      Recombine Surface {s};
    EndFor
  EndFor
EndFor

For j In {0:3}
  For k In {0:2}
    For i In {0:2}
      s = 5000 + i + 3*k + 9*j;
      Transfinite Surface {s};
      Recombine Surface {s};
    EndFor
  EndFor
EndFor

For i In {0:3}
  For k In {0:2}
    For j In {0:2}
      s = 6000 + j + 3*k + 9*i;
      Transfinite Surface {s};
      Recombine Surface {s};
    EndFor
  EndFor
EndFor

// ---------------------------------------------------------------------------
// Volumes
//
// There are 3 x 3 x 3 logical blocks.
// The central block i=1, j=1, k=1 is omitted.
// That omitted block is the solid cube obstacle.
// ---------------------------------------------------------------------------

fluid_vols[] = {};

For k In {0:2}
  For j In {0:2}
    For i In {0:2}

      If (i != 1 || j != 1 || k != 1)

        v  = 9000 + i + 3*j + 9*k;
        sl = 8000 + i + 3*j + 9*k;

        s_z0 = 4000 + i + 3*j + 9*k;
        s_z1 = 4000 + i + 3*j + 9*(k+1);

        s_y0 = 5000 + i + 3*k + 9*j;
        s_y1 = 5000 + i + 3*k + 9*(j+1);

        s_x0 = 6000 + j + 3*k + 9*i;
        s_x1 = 6000 + j + 3*k + 9*(i+1);

        Surface Loop(sl) = {-s_z0, s_z1, s_y0, -s_y1, -s_x0, s_x1};
        Volume(v) = {sl};

        Transfinite Volume {v};

        fluid_vols[] += {v};

      EndIf

    EndFor
  EndFor
EndFor

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

// Cube obstacle faces
obs_ymin() = Surface In BoundingBox {
  x1 - eps, y1 - eps, z1 - eps,
  x2 + eps, y1 + eps, z2 + eps
};

obs_ymax() = Surface In BoundingBox {
  x1 - eps, y2 - eps, z1 - eps,
  x2 + eps, y2 + eps, z2 + eps
};

obs_xmin() = Surface In BoundingBox {
  x1 - eps, y1 - eps, z1 - eps,
  x1 + eps, y2 + eps, z2 + eps
};

obs_xmax() = Surface In BoundingBox {
  x2 - eps, y1 - eps, z1 - eps,
  x2 + eps, y2 + eps, z2 + eps
};

obs_zmin() = Surface In BoundingBox {
  x1 - eps, y1 - eps, z1 - eps,
  x2 + eps, y2 + eps, z1 + eps
};

obs_zmax() = Surface In BoundingBox {
  x1 - eps, y1 - eps, z2 - eps,
  x2 + eps, y2 + eps, z2 + eps
};

Physical Surface("inlet") = {inlet()};
Physical Surface("outlet") = {outlet()};

Physical Surface("wall") = {
  bottom_wall(),
  top_wall(),
  obs_ymin(),
  obs_ymax(),
  obs_xmin(),
  obs_xmax(),
  obs_zmin(),
  obs_zmax()
};

Physical Surface("zmin") = {zmin()};
Physical Surface("zmax") = {zmax()};

Physical Volume("fluid") = {fluid_vols[]};