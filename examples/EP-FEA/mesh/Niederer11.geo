//+
Point(1) = {0, 0, 0,  0.1};
//+
Point(2) = {20, 0, 0, 0.1};
//+
Point(3) = {20, 7, 0, 0.1};
//+
Point(4) = {0, 7, 0,  0.1};
//+
Point(5) = {0, 0, 3,  0.1};
//+
Point(6) = {20, 0, 3, 0.1};
//+
Point(7) = {20, 7, 3, 0.1};
//+
Point(8) = {0, 7, 3,  0.1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 6};
//+
Line(3) = {6, 7};
//+
Line(4) = {7, 3};
//+
Line(5) = {3, 2};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 4};
//+
Line(8) = {4, 1};
//+
Line(9) = {1, 5};
//+
Line(10) = {5, 8};
//+
Line(11) = {5, 6};
//+
Line(12) = {3, 4};
//+
Curve Loop(1) = {6, -10, 11, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 5, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {4, 12, -7, -6};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 8, 1, -5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {2, -11, -9, 1};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {10, 7, 8, 9};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 3, 2, 4, 6, 5};
//+
Volume(1) = {1};
