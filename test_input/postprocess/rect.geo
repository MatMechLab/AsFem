//+
SetFactory("OpenCASCADE");

dx=0.1;

w=0.5;
h=2.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  h,0.0,dx};
Point(4)={0.0,  h,0.0,dx};

//+
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
Physical Surface("block", 9) = {1};
