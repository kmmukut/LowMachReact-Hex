// Unit-cube lid-driven cavity. Generate with:
//   gmsh -3 cavity.geo -o cavity.msh

SetFactory("Built-in");

Lx = 1.0;
Ly = 1.0;
Lz = 1.0/16.0;

nx = 16;
ny = 16;
nz = 1;

Point(1) = {0,  0,  0,  1.0};
Point(2) = {Lx, 0,  0,  1.0};
Point(3) = {Lx, Ly, 0,  1.0};
Point(4) = {0,  Ly, 0,  1.0};
Point(5) = {0,  0,  Lz, 1.0};
Point(6) = {Lx, 0,  Lz, 1.0};
Point(7) = {Lx, Ly, Lz, 1.0};
Point(8) = {0,  Ly, Lz, 1.0};

Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,1};
Line(5)  = {5,6};
Line(6)  = {6,7};
Line(7)  = {7,8};
Line(8)  = {8,5};
Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1}; // zmin
Curve Loop(2) = {5,6,7,8};
Plane Surface(2) = {2}; // zmax
Curve Loop(3) = {1,10,-5,-9};
Plane Surface(3) = {3}; // ymin
Curve Loop(4) = {3,12,-7,-11};
Plane Surface(4) = {4}; // ymax moving lid
Curve Loop(5) = {4,9,-8,-12};
Plane Surface(5) = {5}; // xmin
Curve Loop(6) = {2,11,-6,-10};
Plane Surface(6) = {6}; // xmax

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Transfinite Curve {1,3,5,7} = nx + 1;
Transfinite Curve {2,4,6,8} = ny + 1;
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

