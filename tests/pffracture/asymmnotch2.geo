//+
SetFactory("OpenCASCADE");

dx=1.2;
w=10.0;
dh=1.0e-5;

w=20.0;w1=8.0;
h1=20+1.0;
h2=20-1.0;
dh=1.0e-6;
dh1=3.0;

n1=21;//short notch edge
n2=81;//long in plane edge
n3=5;// short side edges
n4=9;//long side edges

Point(1)={0.0,0.0,0.0,dx};
Point(2)={w*2,0.0,0.0,dx};
Point(3)={w*2,w*2,0.0,dx};
Point(4)={0.0,w*2,0.0,dx};

Point(5)={w*2,h2,0.0,dx};
Point(6)={w*2,h2-dh,0.0,dx};
Point(7)={w*2,h2+dh,0.0,dx};
Point(8)={w*2-8,h2,0.0,dx};

Point( 9)={0.0,h1,0.0,dx};
Point(10)={0.0,h1-dh,0.0,dx};
Point(11)={0.0,h1+dh,0.0,dx};
Point(12)={8.0,h1,0.0,dx};

Point(13)={0.0,h1+dh1,0.0,dx};
Point(14)={w*2,h1+dh1,0.0,dx};

Point(15)={0.0,h2-dh1,0.0,dx};
Point(16)={w*2,h2-dh1,0.0,dx};//+
Line(1) = {1, 2};
//+
Line(2) = {2, 16};
//+
Line(3) = {16, 6};
//+
Line(4) = {6, 8};
//+
Line(5) = {8, 7};
//+
Line(6) = {7, 14};
//+
Line(7) = {14, 3};
//+
Line(8) = {3, 4};
//+
Line(9) = {4, 13};
//+
Line(10) = {13, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 10};
//+
Line(13) = {10, 15};
//+
Line(14) = {15, 1};
//+
Line(15) = {13, 14};
//+
Line(16) = {15, 16};
//+
Curve Loop(1) = {14, 1, 2, -16};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, 16, 3, 4, 5, 6, -15, 10, 11, 12};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 9, 15, 7};
//+
Plane Surface(3) = {3};
//+
Physical Curve("left") = {9, 10, 13, 14};
//+
Physical Curve("bottom") = {1};
//+
Physical Curve("right") = {2, 3, 6, 7};
//+
Physical Curve("top") = {8};
//+
Physical Surface("block") = {1, 2, 3};
//+
Transfinite Curve {11, 12, 5, 4} = n1 Using Progression 1;
//+
Transfinite Curve {15, 16} = n2 Using Progression 1;
//+
Transfinite Curve {10, 3} = n3 Using Progression 1;
//+
Transfinite Curve {13, 6} = n4 Using Progression 1;
