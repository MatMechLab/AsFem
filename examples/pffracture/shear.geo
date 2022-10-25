//+
SetFactory("OpenCASCADE");

w=1.0;

dx=0.025;
dx1=0.0015;
dx2=0.0008;

dw=1.0e-5;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx1};
Point(3)={  w,  w,0.0,dx};
Point(4)={0.0,  w,0.0,dx};

Point(5)={w/2,w/2,0.0,dx2};

Point(6)={0.0,w/2-dw,0.0,dx};
Point(7)={0.0,w/2+dw,0.0,dx};

//+
Point(8) = {0.85, 0.1, 0, dx1};
//+
Point(9) = {0.75, 0.2, 0, dx1};
//+
Point(10) = {0.55, 0.4, 0, dx1};
//+
Point(11) = {0.65, 0.3, 0, dx1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 7};
//+
Line(5) = {7, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 1};
//+
Spline(8) = {5, 10, 11, 9, 8, 2};
//+
Curve Loop(1) = {7, 1, -8, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 5, 8, 2, 3};
//+
Plane Surface(2) = {2};
//+
Physical Curve("left", 9) = {4, 7};
//+
Physical Curve("bottom", 10) = {1};
//+
Physical Curve("right", 11) = {2};
//+
Physical Curve("top", 12) = {3};
//+
Physical Surface("block", 13) = {2, 1};
