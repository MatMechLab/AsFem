//+
SetFactory("OpenCASCADE");

dx=0.02;
dw=1.0e-6;
dw1=0.012;

n1=211;
n2=8;
n3=16;

dx=0.04;
dw1=0.0095;
n1=501;
n2=11;
n3=21;
n4=81;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={1.0,0.0,0.0,dx};
Point(3)={1.0,1.0,0.0,dx};
Point(4)={0.0,1.0,0.0,dx};

Point(5)={0.0,0.5-dw,0.0,dx};
Point(6)={0.0,0.5+dw,0.0,dx};

Point(7)={0.5,    0.5,0.0,dx/50.0};
Point(8)={0.5,0.5-dw1,0.0,dx/50.0};
Point(9)={0.5,0.5+dw1,0.0,dx/50.0};

Point(10)={1.0,0.5-dw1,0.0,dx};
Point(11)={1.0,0.5+dw1,0.0,dx};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 10};
//+
Line(3) = {10, 8};
//+
Line(4) = {8, 7};
//+
Line(5) = {7, 9};
//+
Line(6) = {9, 11};
//+
Line(7) = {10, 11};
//+
Line(8) = {11, 3};
//+
Line(9) = {3, 4};
//+
Line(10) = {4, 6};
//+
Line(11) = {6, 7};
//+
Line(12) = {7, 5};
//+
Line(13) = {5, 1};
//+
Curve Loop(1) = {13, 1, 2, 3, 4, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 5, 6, -7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 9, 10, 11, 5, 6};
//+
Plane Surface(3) = {3};
//+
Physical Curve("left", 14) = {10, 13};
//+
Physical Curve("bottom", 15) = {1};
//+
Physical Curve("right", 16) = {2, 7, 8};
//+
Physical Curve("top", 17) = {9};
//+
Physical Surface("block", 18) = {1, 2, 3};
//+
Transfinite Curve {6, 3} = n1 Using Progression 1;
//+
//Transfinite Curve {5, 4} = n2 Using Progression 1;
//+
Transfinite Curve {7} = n3 Using Progression 1;//+
Transfinite Curve {9} = n4 Using Progression 1;
