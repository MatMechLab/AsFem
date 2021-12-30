//+
SetFactory("OpenCASCADE");

dx=0.01;
dw=1.0e-6;

n=201;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={1.0,0.0,0.0,dx};
Point(3)={0.5,0.5,0.0,dx};
Point(4)={1.0,1.0,0.0,dx};
Point(5)={0.0,1.0,0.0,dx};

Point(6)={0.0,0.5+dw,0.0,dx};
Point(7)={0.0,0.5-dw,0.0,dx};//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 3};
//+
Line(6) = {3, 7};
//+
Line(7) = {7, 1};
//+
Line(8) = {3, 2};
//+
Curve Loop(1) = {7, 1, -8, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 5, 8, 2};
//+
Plane Surface(2) = {2};
//+
Physical Curve("bottom") = {1};
//+
Physical Curve("right") = {2};
//+
Physical Curve("top") = {3};
//+
Physical Curve("left") = {4, 7};
//+
//Transfinite Curve {8} = n Using Progression 1;
//+
Physical Surface("block") = {1, 2};