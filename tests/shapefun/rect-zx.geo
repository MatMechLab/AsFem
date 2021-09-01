//+
SetFactory("OpenCASCADE");
dx=0.1;

w=1.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,0.0,w,dx};
Point(4)={0.0,0.0,w,dx};//+
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {1};
//+
Physical Line("right") = {3};
//+
Physical Surface("block") = {1};
