//+
SetFactory("OpenCASCADE");
dx=0.1;
y=1.0;z=2.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={0.0,  y,0.0,dx};
Point(3)={0.0,  y,  z,dx};
Point(4)={0.0,0.0,  z,dx};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("bottom") = {4};
//+
Physical Line("top") = {2};
//+
Physical Line("back") = {1};
//+
Physical Line("front") = {3};
//+
Physical Surface("block") = {1};
