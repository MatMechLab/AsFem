//+
SetFactory("OpenCASCADE");

dx=0.5;

w=1.0;
Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  w,0.0,dx};
Point(4)={0.0,  w,0.0,dx};

Point(5)={0.0,0.0,  w,dx};
Point(6)={  w,0.0,  w,dx};
Point(7)={  w,  w,  w,dx};
Point(8)={0.0,  w,  w,dx};//+
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
Curve Loop(1) = {4, -9, -8, 12};
//+
Plane Surface(1) = {1};

//+
Curve Loop(2) = {10, 2, -11, -6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {9, 1, -10, -5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, -3, -11, 7};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {8, 5, 6, 7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {4, 1, 2, 3};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {5, 1, 6, 3, 2, 4};
//+
Volume(1) = {1};
//+
Physical Surface("left", 13) = {1};
//+
Physical Surface("right", 14) = {2};
//+
Physical Surface("bottom", 15) = {3};
//+
Physical Surface("top", 16) = {4};
//+
Physical Surface("back", 17) = {6};
//+
Physical Surface("front", 18) = {5};
//+
Physical Volume("matrix", 19) = {1};
