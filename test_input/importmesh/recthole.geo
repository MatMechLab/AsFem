//+
SetFactory("OpenCASCADE");

dx=0.2;
w=1.0;
r=0.2;

n=26;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  w,0.0,dx};
Point(4)={0.0,  w,0.0,dx};

Point(5)={w/2,w/2+0.1,0.0,dx};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {0.5, 0.5, 0, r, 0, 2*Pi};
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
Physical Curve("left", 10) = {6};
//+
Physical Curve("bottom", 11) = {8};
//+
Physical Curve("right", 12) = {9};
//+
Physical Curve("top", 13) = {7};
//+
Physical Surface("matrix", 14) = {1};
//+
Physical Surface("inclusion", 15) = {2};
//+
Transfinite Curve {5} = n Using Progression 1;
