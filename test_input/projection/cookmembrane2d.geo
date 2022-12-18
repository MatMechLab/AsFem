//+
SetFactory("OpenCASCADE");

dx=2.5;
w=48.0;
h1=44.0;
h2=16.0;

Point(1)={0.0,  0.0,0.0,dx};
Point(2)={  w,   h1,0.0,dx};
Point(3)={  w,h1+h2,0.0,dx};
Point(4)={0.0,   h1,0.0,dx/2.0};//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("left", 5) = {4};
//+
Physical Curve("bottom", 6) = {1};
//+
Physical Curve("right", 7) = {2};
//+
Physical Curve("top", 8) = {3};
//+
Physical Surface("matrix", 9) = {1};
//+
Physical Point("A", 10) = {3};
