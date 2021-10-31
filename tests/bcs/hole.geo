//+
SetFactory("OpenCASCADE");

dx=0.02;
r=0.4;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={1.0,0.0,0.0,dx};
Point(3)={1.0,1.0,0.0,dx};
Point(4)={0.0,1.0,0.0,dx};//+

Point(5)={  r,0.0,0.0,dx};
Point(6)={0.0,  r,0.0,dx};
//+
Circle(1) = {6, 1, 5};
//+
Line(2) = {5, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 6};
//+
Curve Loop(1) = {5, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("surface", 6) = {1};
//+
Physical Surface("block", 7) = {1};
