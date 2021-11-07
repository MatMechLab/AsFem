//+
SetFactory("OpenCASCADE");

dx=0.006;

r=0.085;

Mesh.MinimumCircleNodes=101;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={1.0,0.0,0.0,dx};
Point(3)={1.0,1.0,0.0,dx};
Point(4)={0.0,1.0,0.0,dx};//+
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
Circle(5) = {0.16, 0.72, 0, r, 0, 2*Pi};
//+
Circle(6) = {0.44, 0.6, 0, r, 0, 2*Pi};
//+
Circle(7) = {0.7, 0.7, 0, r, 0, 2*Pi};
//+
Circle(8) = {0.85, 0.42, 0, r, 0, 2*Pi};
//+
Circle(9) = {0.62, 0.2, -0, r, 0, 2*Pi};
//+
Circle(10) = {0.3, 0.3, -0, r, 0, 2*Pi};
//+
Circle(11) = {0.15, 0.45, -0, r, 0, 2*Pi};
//+
Circle(12) = {0.4, 0.83, 0, r, 0, 2*Pi};
//+
Circle(13) = {0.6, 0.4, 0, r, 0, 2*Pi};
//+
Circle(14) = {0.84, 0.19, 0, r, 0, 2*Pi};
//+
Circle(15) = {0.16, 0.15, 0, r, 0, 2*Pi};
//+
Circle(16) = {0.84, 0.85, 0, r, 0, 2*Pi};
//+
Curve Loop(2) = {5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {16};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {6};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {11};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {10};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {13};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {8};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {9};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {14};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {15};
//+
Plane Surface(13) = {13};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{13}; Surface{9}; Surface{10}; Surface{12}; Surface{11}; Delete; }
//+
Physical Curve("left", 21) = {17};
//+
Physical Curve("bottom", 22) = {19};
//+
Physical Curve("right", 23) = {20};
//+
Physical Curve("top", 24) = {18};
//+
Physical Surface("block", 25) = {1};
