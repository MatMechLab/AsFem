//+
SetFactory("OpenCASCADE");

dx0=0.025;
dx1=0.0035;

w=1.0;
h=1.0;


// coarse
dh=2.0e-4;
dx0=0.1;
dx1=0.01;

dx0=0.025;
dx1=0.004;

//dh=1.0e-4;
//dx0=0.075;
//dx1=0.006;

// for fine mesh
//dh=1.0e-4;
//dx0=0.025;
//dx1=0.0025;

dw=0.02;
dz=0.05;
n1=61;

Mesh.AngleToleranceFacetOverlap=0.01;

Point(1)={0.0,0.0,0.0,dx0};
Point(2)={  w,0.0,0.0,dx0};
Point(3)={  w,  h,0.0,dx0};
Point(4)={0.0,  h,0.0,dx0};
Point(5)={0.0,h/2+dh,0.0,dx0};
Point(6)={0.0,h/2-dh,0.0,dx0};
Point(7)={w/2,h/2   ,0.0,dx1/2};
// for refine region
Point(8)={w/2,h/2+dw,0.0,dx1};
Point(9)={w/2,h/2-dw,0.0,dx1};
Point(10)={w,h/2+dw,0.0,dx1};
Point(11)={w,h/2-dw,0.0,dx1};
// for top front layer
Point(12)={0.0,0.0,dz,dx0};
Point(13)={  w,0.0,dz,dx0};
Point(14)={  w,  h,dz,dx0};
Point(15)={0.0,  h,dz,dx0};
Point(16)={0.0,h/2+dh,dz,dx0};
Point(17)={0.0,h/2-dh,dz,dx0};
Point(18)={w/2,h/2   ,dz,dx1/2};
// for refine region
Point(19)={w/2,h/2+dw,dz,dx1};
Point(20)={w/2,h/2-dw,dz,dx1};
Point(21)={w,h/2+dw,dz,dx1};
Point(22)={w,h/2-dw,dz,dx1};//+
Line(1) = {12, 13};
//+
Line(2) = {13, 22};
//+
Line(3) = {22, 21};
//+
Line(4) = {21, 14};
//+
Line(5) = {14, 15};
//+
Line(6) = {15, 16};
//+
Line(7) = {16, 18};
//+
Line(8) = {18, 17};
//+
Line(9) = {17, 12};
//+
Line(10) = {18, 20};
//+
Line(11) = {20, 22};
//+
Line(12) = {18, 19};
//+
Line(13) = {19, 21};
//+
Line(14) = {1, 2};
//+
Line(15) = {2, 11};
//+
Line(16) = {11, 10};
//+
Line(17) = {10, 3};
//+
Line(18) = {3, 4};
//+
Line(19) = {4, 5};
//+
Line(20) = {5, 7};
//+
Line(21) = {7, 6};
//+
Line(22) = {6, 1};
//+
Line(23) = {7, 8};
//+
Line(24) = {8, 10};
//+
Line(25) = {7, 9};
//+
Line(26) = {9, 11};
//+
Line(27) = {1, 12};
//+
Line(28) = {6, 17};
//+
Line(29) = {5, 16};
//+
Line(30) = {7, 18};
//+
Line(31) = {9, 20};
//+
Line(32) = {8, 19};
//+
Line(33) = {4, 15};
//+
Line(34) = {3, 14};
//+
Line(35) = {10, 21};
//+
Line(36) = {11, 22};
//+
Line(37) = {2, 13};
//+
Curve Loop(1) = {27, 1, -37, -14};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, -27, -22, 28};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, -28, -21, 30};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, -30, -20, 29};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {11, -36, -26, 31};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {15, 36, -2, -37};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {14, 15, -26, -25, 21, 22};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {9, 1, 2, -11, -10, 8};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {31, -10, -30, 25};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {12, -32, -23, 30};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {6, -29, -19, 33};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {13, -35, -24, 32};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {13, -3, -11, -10, 12};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {16, -24, -23, 25, 26};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {3, -35, -16, 36};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {7, 12, 13, 4, 5, 6};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {17, 18, 19, 20, 23, 24};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {17, 34, -4, -35};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {5, -33, -18, 34};
//+
Plane Surface(19) = {19};
//+
Surface Loop(1) = {8, 2, 1, 6, 7, 3, 9, 5};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {16, 4, 17, 18, 19, 11, 12, 10};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {13, 15, 14, 5, 9, 10, 12};
//+
Volume(3) = {3};
//+
Physical Surface("left", 38) = {11, 2};
//+
Physical Surface("bottom", 39) = {1};
//+
Physical Surface("right", 40) = {6, 15, 18};
//+
Physical Surface("top", 41) = {19};
//+
Physical Surface("back", 42) = {17, 14, 7};
//+
Physical Surface("front", 43) = {16, 13, 8};
//+
Physical Volume("block", 44) = {2, 3, 1};
//+
Transfinite Curve {30} = n1 Using Progression 1;
