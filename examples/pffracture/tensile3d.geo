//+
SetFactory("OpenCASCADE");

dx=0.0225;
dx1=0.0035;

w=1.0;
w1=0.5;
h=0.5;
z=0.025;
dh=6.0e-4;

n1=601;// for long edge
n2=81;// for short edge

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  h,0.0,dx};
Point(4)={0.0,  h,0.0,dx};
Point(5)={0.0,h/2+dh,0.0,dx};
Point(6)={0.0,h/2-dh,0.0,dx};
Point(7)={ w1,h/2   ,0.0,dx1};
Point(8)={  w,h/2   ,0.0,dx1};

// front layer

Point(9) ={0.0,0.0,z,dx};
Point(10)={  w,0.0,z,dx};
Point(11)={  w,  h,z,dx};
Point(12)={0.0,  h,z,dx};
Point(13)={0.0,h/2+dh,z,dx};
Point(14)={0.0,h/2-dh,z,dx};
Point(15)={ w1,h/2   ,z,dx1};
Point(16)={  w,h/2   ,z,dx1};//+
Line(1) = {9, 10};
//+
Line(2) = {10, 16};
//+
Line(3) = {16, 11};
//+
Line(4) = {11, 12};
//+
Line(5) = {12, 13};
//+
Line(6) = {13, 15};
//+
Line(7) = {15, 14};
//+
Line(8) = {14, 9};
//+
Line(9) = {15, 16};
//+
Line(10) = {1, 2};
//+
Line(11) = {2, 8};
//+
Line(12) = {8, 3};
//+
Line(13) = {3, 4};
//+
Line(14) = {4, 5};
//+
Line(15) = {5, 7};
//+
Line(16) = {7, 6};
//+
Line(17) = {6, 1};
//+
Line(18) = {7, 8};
//+
Line(19) = {1, 9};
//+
Line(20) = {2, 10};
//+
Line(21) = {8, 16};
//+
Line(22) = {3, 11};
//+
Line(23) = {4, 12};
//+
Line(24) = {5, 13};
//+
Line(25) = {7, 15};
//+
Line(26) = {6, 14};
//+
Curve Loop(1) = {10, 20, -1, -19};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {17, 19, -8, -26};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 2, -9, 7, 8};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11, 21, -2, -20};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {18, -11, -10, -17, -16};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {16, 26, -7, -25};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {15, 25, -6, -24};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {18, 21, -9, -25};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {12, 22, -3, -21};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {14, 24, -5, -23};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {13, 23, -4, -22};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {9, 3, 4, 5, 6};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {13, 14, 15, 18, 12};
//+
Plane Surface(13) = {13};
//+
Surface Loop(1) = {1, 5, 4, 3, 6, 2, 8};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {13, 11, 10, 7, 12, 9, 8};
//+
Volume(2) = {2};
//+
Physical Surface("left", 27) = {10, 2};
//+
Physical Surface("right", 28) = {4, 9};
//+
Physical Surface("bottom", 29) = {1};
//+
Physical Surface("top", 30) = {11};
//+
Physical Surface("back", 31) = {5, 13};
//+
Physical Surface("front", 32) = {3, 12};
//+
Physical Volume("matrix", 33) = {1, 2};
//+
//Transfinite Curve {18, 9} = n1 Using Progression 1;
//+
//Transfinite Curve {25, 21} = n2 Using Progression 1;
