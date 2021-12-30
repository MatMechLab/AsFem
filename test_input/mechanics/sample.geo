SetFactory("OpenCASCADE");
dx=1.0;//+

w1=40.0;w2=100.0;
r=12.0;
r1=20.0-r;

Point(1)={0.0,r1,0.0,dx};

Point(2)={w2/2,r1,0.0,dx};
Point(3)={w2/2,20,0.0,dx};
Point(4)={w2/2+r,20,0.0,dx};
Point(5)={w2/2+r+w1,20,0.0,dx};


Point(6)={-w2/2,r1,0.0,dx};
Point(7)={-w2/2,20,0.0,dx};
Point(8)={-(w2/2+r),20,0.0,dx};
Point(9)={-(w2/2+r+w1),20,0.0,dx};


Point(10)={w2/2,-r1,0.0,dx};
Point(11)={w2/2,-20,0.0,dx};
Point(12)={w2/2+r,-20,0.0,dx};
Point(13)={w2/2+r+w1,-20,0.0,dx};


Point(14)={-w2/2,-r1,0.0,dx};
Point(15)={-w2/2,-20,0.0,dx};
Point(16)={-(w2/2+r),-20,0.0,dx};
Point(17)={-(w2/2+r+w1),-20,0.0,dx};//+
Line(1) = {9, 8};
//+
Line(2) = {6, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 13};
//+
Line(6) = {13, 12};
//+
Line(7) = {10, 14};
//+
Line(8) = {16, 17};
//+
Line(9) = {17, 9};
//+
Circle(10) = {6, 7, 8};
//+
Circle(11) = {16, 15, 14};
//+
Circle(12) = {4, 3, 2};
//+
Circle(13) = {10, 11, 12};
//+
Curve Loop(1) = {9, 1, -10, 2, 3, -12, 4, 5, 6, -13, 7, -11, 8};
//+
Plane Surface(1) = {1};
//+
Physical Curve("left") = {9};
//+
Physical Curve("right") = {5};
//+
Physical Surface("block") = {1};
