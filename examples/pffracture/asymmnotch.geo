//+
SetFactory("OpenCASCADE");

dx=0.05;
w=10.0;
dh=1.0e-5;

n1=151;//top-two lines
n2=81;// bottom circle
n3=51;// left-right lines
n4=61;// bottom-notch two lines

Point(1)={0.0,0.0,0.0,dx};
Point(2)={1.0,0.0,0.0,dx};
Point(3)={4.0-dh,0.0,0.0,dx};
Point(4)={4.0+dh,0.0,0.0,dx};
Point(5)={w+w-1.0,0.0,0.0,dx};
Point(6)={w+w,0.0,0.0,dx};
Point(7)={w*2,0.0,0.0,dx};
Point(8)={w*2,8.0,0.0,0.0,dx};
Point(9)={0.0,8.0,0.0,dx};
Point(10)={4.0,1.0,0.0,0.0,dx};
Point(11)={w,8.0,0.0,0.0,dx};

//+
Circle(1) = {6.0,2.75,0.0,0.25,0, 2*Pi};
//+
Circle(2) = {6.0,2.75+2,0.0,0.25,0, 2*Pi};
//+
Circle(3) = {6.0,2.75+2+2,0.0,0.25,0, 2*Pi};//+
//+
Line(4) = {1, 2};
//+
Line(5) = {2, 3};
//+
Line(6) = {3, 10};
//+
Line(7) = {10, 4};
//+
Line(8) = {4, 5};
//+
Line(9) = {5, 6};
//+
Line(10) = {6, 8};
//+
Line(11) = {8, 11};
//+
Line(12) = {11, 9};
//+
Line(13) = {9, 1};
//+
Curve Loop(1) = {13, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1};
//+
Plane Surface(4) = {4};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Delete; }
//+
Physical Point("leftpoint") = {18};
//+
Physical Point("rightpoint") = {23};
//+
Physical Point("toppoint") = {17};
//+
Physical Curve("top") = {7, 5};
//+
Physical Surface("block") = {1};
//+
Transfinite Curve {5, 7} = n1 Using Progression 1;
//+
Transfinite Curve {1,2} = n2 Using Progression 1;
//+
Transfinite Curve {3} = n2 Using Progression 1;
//+
Transfinite Curve {4} = 81 Using Progression 1;
Transfinite Curve {9} = n3 Using Progression 1;
//+
Transfinite Curve {10, 12} = n4 Using Progression 1;
