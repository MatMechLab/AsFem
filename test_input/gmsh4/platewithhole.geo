//+
SetFactory("OpenCASCADE");

dx=0.08;
r=1.0;
w=5.0;

n=31;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  w,0.0,dx};
Point(4)={0.0,  w,0.0,dx};

Point(5)={  r,0.0,0.0,dx};
Point(6)={0.0,  r,0.0,dx};//+
Line(1) = {5, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 6};
//+
Circle(5) = {6, 1, 5};
//+
Curve Loop(1) = {4, 5, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Point("leftpoint") = {6};
//+
Physical Point("bottompoint") = {5};
//+
Physical Curve("left") = {4};
//+
Physical Curve("bottom") = {1};
//+
Physical Curve("right") = {2};
//+
Physical Curve("top") = {3};
//+
Physical Surface("block") = {1};
//+
Transfinite Curve {5} = n Using Progression 1;
