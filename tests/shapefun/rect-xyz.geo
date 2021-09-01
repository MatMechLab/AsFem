//+
SetFactory("OpenCASCADE");
dx=0.2;

w=1.0;

Point(1)={w,w,0,dx};
Point(2)={w,0,w,dx};
Point(3)={0,w,w,dx};
Point(4)={0,0,2*w,dx};
//+
Line(1) = {4, 2};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 4};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {4};
//+
Physical Line("right") = {2};
//+
Physical Surface("block") = {1};
