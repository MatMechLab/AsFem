//+
SetFactory("OpenCASCADE");

dx=0.4;
w=1.0;
r=0.2;

n=6;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  w,0.0,dx};
Point(4)={0.0,  w,0.0,dx};

Point(5)={w/2,w/2,0.0,dx};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {w/2, w/2, 0, r, 0, 2*Pi};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5};
//+
Plane Surface(2) = {2};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Curve Loop(1) = {5};
//+
Plane Surface(2) = {1};
//+
Physical Curve("left", 6) = {1};
//+
Physical Curve("bottom", 7) = {3};
//+
Physical Curve("right", 8) = {4};
//+
Physical Curve("top", 9) = {2};
//+
Physical Surface("matrix", 10) = {1};
//+
Physical Surface("inclusion", 11) = {2};
//+
Transfinite Curve {5} = n Using Progression 1;
