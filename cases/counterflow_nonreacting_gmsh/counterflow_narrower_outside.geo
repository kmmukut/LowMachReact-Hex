// Thin 2D counterflow slab for validation against Cantera CounterflowDiffusionFlame:
//   - Streamwise x between opposed inlets (xmin / xmax).
//   - Resolved transverse y direction with open outlets (ymin / ymax).
//   - Exactly ONE hex cell in z (nz=1), with zmin/zmax as symmetry planes.
//   - Center-clustered structured hex mesh near x=Lx/2, y=Ly/2.
//
// Generate native mesh:
//   ./make_mesh.sh
//
// Physical patch names: xmin, xmax, ymin, ymax, zmin, zmax.

SetFactory("Built-in");

Lx = 2.0;
Ly = 1.0;
Lz = 0.02;

// Total cell counts.
// Keep these modest while testing runtime.
nx = 64;
ny = 32;
nz = 1;

// Center clustering strength.
// Smaller values generally give stronger clustering toward the curve center.
// Good starting values:
//   0.20 to 0.35 for x
//   0.25 to 0.45 for y
xBump = 0.25;
yBump = 0.35;

Point(1) = {0,  0,  0,  1.0};
Point(2) = {Lx, 0,  0,  1.0};
Point(3) = {Lx, Ly, 0,  1.0};
Point(4) = {0,  Ly, 0,  1.0};
Point(5) = {0,  0,  Lz, 1.0};
Point(6) = {Lx, 0,  Lz, 1.0};
Point(7) = {Lx, Ly, Lz, 1.0};
Point(8) = {0,  Ly, Lz, 1.0};

Line(1)  = {1,2};  // x direction, y=0, z=0
Line(2)  = {2,3};  // y direction, x=Lx, z=0
Line(3)  = {3,4};  // x direction, y=Ly, z=0
Line(4)  = {4,1};  // y direction, x=0, z=0
Line(5)  = {5,6};  // x direction, y=0, z=Lz
Line(6)  = {6,7};  // y direction, x=Lx, z=Lz
Line(7)  = {7,8};  // x direction, y=Ly, z=Lz
Line(8)  = {8,5};  // y direction, x=0, z=Lz
Line(9)  = {1,5};  // z direction
Line(10) = {2,6};  // z direction
Line(11) = {3,7};  // z direction
Line(12) = {4,8};  // z direction

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1}; // zmin

Curve Loop(2) = {5,6,7,8};
Plane Surface(2) = {2}; // zmax

Curve Loop(3) = {1,10,-5,-9};
Plane Surface(3) = {3}; // ymin

Curve Loop(4) = {3,12,-7,-11};
Plane Surface(4) = {4}; // ymax

Curve Loop(5) = {4,9,-8,-12};
Plane Surface(5) = {5}; // xmin

Curve Loop(6) = {2,11,-6,-10};
Plane Surface(6) = {6}; // xmax

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

// Center-clustered transfinite spacing.
// x-direction curves: refine near x = Lx/2.
Transfinite Curve {1,3,5,7} = nx + 1 Using Bump xBump;

// y-direction curves: refine near y = Ly/2.
// This helps the extracted centerline and the core opposed-flow region.
Transfinite Curve {2,4,6,8} = ny + 1 Using Bump yBump;

// z remains one cell.
Transfinite Curve {9,10,11,12} = nz + 1;

Transfinite Surface "*";
Transfinite Volume {1};
Recombine Surface "*";

Physical Surface("zmin") = {1};
Physical Surface("zmax") = {2};
Physical Surface("ymin") = {3};
Physical Surface("ymax") = {4};
Physical Surface("xmin") = {5};
Physical Surface("xmax") = {6};
Physical Volume("fluid") = {1};