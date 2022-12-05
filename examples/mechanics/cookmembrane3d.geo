//+
SetFactory("OpenCASCADE");

dx=4.0;
w=48.0;
h1=44.0;
h2=16.0;

z=1.0;

Point(1)={0.0,  0.0,z/2.0,dx};
Point(2)={  w,   h1,z/2.0,dx};
Point(3)={  w,h1+h2,z/2.0,dx};
Point(4)={0.0,   h1,z/2.0,dx/2};//+

Point(5)={0.0,  0.0,-z/2.0,dx};
Point(6)={  w,   h1,-z/2.0,dx};
Point(7)={  w,h1+h2,-z/2.0,dx};
Point(8)={0.0,   h1,-z/2.0,dx/2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {5, 1};
//+
Line(10) = {6, 2};
//+
Line(11) = {7, 3};
//+
Line(12) = {8, 4};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {4, -9, -8, 12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, -11, -6, 10};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {3, -12, -7, 11};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {1, -10, -5, 9};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {2, 6, 1, 3, 5, 4};
//+
Volume(1) = {1};
//+
Physical Surface("left", 13) = {3};
//+
Physical Surface("right", 14) = {4};
//+
Physical Surface("bottom", 15) = {6};
//+
Physical Surface("top", 16) = {5};
//+
Physical Surface("back", 17) = {2};
//+
Physical Surface("front", 18) = {1};
//+
Physical Volume("matrix", 19) = {1};
