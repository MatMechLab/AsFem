//+
SetFactory("OpenCASCADE");

dx=0.05;
w=4.0;
dw=1.0e-1;
dw1=0.06;
dw2=0.05;

n1=261; // for middle line
n2=301; // for middle-two-side lines
n3=6;  // for top-two-lines
n4=3;   // for bottom-two-lines
n5=41;  // for bottom-two-slope line




Point(1)={0.0,0.0,0.0,dx};
Point(2)={w-dw-dw1,0.0,0.0,dx};
Point(3)={w-dw    ,0.0,0.0,dx};
Point(4)={w       ,0.4,0.0,dx};
Point(5)={w+dw    ,0.0,0.0,dx};
Point(6)={w+dw+dw1,0.0,0.0,dx};

Point(7)={w+w,0.0,0.0,dx};
Point(8)={w+w,2.0,0.0,dx};
Point(9)={w  ,2.0,0.0,dx};
Point(10)={0.0,2.0,0.0,dx};

Point(11)={w-dw2,2.0,0.0,dx};
Point(12)={w+dw2,2.0,0.0,dx};//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 12};
//+
Line(9) = {12, 9};
//+
Line(10) = {9, 11};
//+
Line(11) = {11, 10};
//+
Line(12) = {10, 1};
//+
Line(13) = {2, 11};
//+
Line(14) = {6, 12};
//+
Line(15) = {4, 9};
//+
Curve Loop(1) = {12, 1, 13, 11};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, -10, -15, -3, -2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {15, -9, -14, -5, -4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, -14, 6, 7};
//+
Plane Surface(4) = {4};
//+
Physical Point("leftpoint") = {1};
//+
Physical Point("rightpoint") = {7};
//+
Physical Point("toppoint") = {9};
//+
Physical Curve("left") = {12};
//+
Physical Curve("right") = {7};
//+
Physical Curve("top") = {11, 8};
//+
Physical Curve("topload") = {10, 9};
//+
Physical Surface("block") = {1, 2, 3, 4};
//+
Transfinite Curve {15} = n1 Using Progression 1;// for middle line
//+
Transfinite Curve {13, 14} = n2 Using Progression 1;// for middle-two-side lines
//+
Transfinite Curve {10, 9} = n3 Using Progression 1;// for top-two-lines
//+
Transfinite Curve {2, 5} = n4 Using Progression 1;// for bottom-two-lines
//+
Transfinite Curve {3, 4} = n5 Using Progression 1;// for bottom-two-slope line
