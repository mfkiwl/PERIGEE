// Gmsh project created on Tue Dec  8 20:27:17 2020
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 1, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {0, 0, 1, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 1};
Line(4) = {3, 2};
Line(5) = {2, 4};
Line(6) = {2, 1};
//+
Curve Loop(1) = {6, -3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 2, 3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, 2, 4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1, -5, 6};
//+
Plane Surface(4) = {4};
//+
Surface Loop(1) = {3, 4, 2, 1};
//+
Volume(1) = {1};
//+
Physical Surface("back") = {4};
//+
Physical Surface("left") = {3};
//+
Physical Surface("diagface") = {2};
//+
Physical Surface("bottom") = {1};
