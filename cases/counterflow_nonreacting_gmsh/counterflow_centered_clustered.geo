// Thin 2D counterflow slab with smaller cells near center.
// Axis-aligned cuboid hex cells, nz = 1.

SetFactory("Built-in");

Lx = 2.0;
Ly = 1.0;
Lz = 0.02;

nx = 64;
ny = 32;
nz = 1;

// Split at center.
xc = 0.5 * Lx;
yc = 0.5 * Ly;

nx_half = nx / 2;
ny_half = ny / 2;

// Center-clustering ratios.
// Larger values => stronger clustering toward center.
// Try 1.03 mild, 1.06 moderate, 1.10 strong.
rx = 1.1;
ry = 1.1;

eps = 1.0e-8;

// Points on z = 0 plane.
Point(1) = {0,  0,  0,  1.0};
Point(2) = {xc, 0,  0,  1.0};
Point(3) = {Lx, 0,  0,  1.0};

Point(4) = {0,  yc, 0,  1.0};
Point(5) = {xc, yc, 0,  1.0};
Point(6) = {Lx, yc, 0,  1.0};

Point(7) = {0,  Ly, 0,  1.0};
Point(8) = {xc, Ly, 0,  1.0};
Point(9) = {Lx, Ly, 0,  1.0};

// Horizontal x-lines.
Line(1) = {1,2};
Line(2) = {2,3};

Line(3) = {4,5};
Line(4) = {5,6};

Line(5) = {7,8};
Line(6) = {8,9};

// Vertical y-lines.
Line(7)  = {1,4};
Line(8)  = {4,7};

Line(9)  = {2,5};
Line(10) = {5,8};

Line(11) = {3,6};
Line(12) = {6,9};

// Four rectangular surfaces.
Curve Loop(1) = {1, 9, -3, -7};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 11, -4, -9};
Plane Surface(2) = {2};

Curve Loop(3) = {3, 10, -5, -8};
Plane Surface(3) = {3};

Curve Loop(4) = {4, 12, -6, -10};
Plane Surface(4) = {4};

// x-direction spacing.
// Left half: large near x=0, small near x=xc.
Transfinite Curve {1,3,5} = nx_half + 1 Using Progression 1/rx;

// Right half: small near x=xc, large near x=Lx.
Transfinite Curve {2,4,6} = nx_half + 1 Using Progression rx;

// y-direction spacing.
// Bottom half: large near y=0, small near y=yc.
Transfinite Curve {7,9,11} = ny_half + 1 Using Progression 1/ry;

// Top half: small near y=yc, large near y=Ly.
Transfinite Curve {8,10,12} = ny_half + 1 Using Progression ry;

// Structured quad surfaces.
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};

Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};

// Extrude one cell in z and recombine into hexes.
Extrude {0, 0, Lz} {
  Surface{1,2,3,4};
  Layers{nz};
  Recombine;
}

// Physical patch names.
xmin[] = Surface In BoundingBox{-eps, -eps, -eps, eps, Ly+eps, Lz+eps};
xmax[] = Surface In BoundingBox{Lx-eps, -eps, -eps, Lx+eps, Ly+eps, Lz+eps};

ymin[] = Surface In BoundingBox{-eps, -eps, -eps, Lx+eps, eps, Lz+eps};
ymax[] = Surface In BoundingBox{-eps, Ly-eps, -eps, Lx+eps, Ly+eps, Lz+eps};

zmin[] = Surface In BoundingBox{-eps, -eps, -eps, Lx+eps, Ly+eps, eps};
zmax[] = Surface In BoundingBox{-eps, -eps, Lz-eps, Lx+eps, Ly+eps, Lz+eps};

Physical Surface("xmin") = {xmin[]};
Physical Surface("xmax") = {xmax[]};

Physical Surface("ymin") = {ymin[]};
Physical Surface("ymax") = {ymax[]};

Physical Surface("zmin") = {zmin[]};
Physical Surface("zmax") = {zmax[]};

Physical Volume("fluid") = Volume{:};