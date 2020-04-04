dx=0.008;

w=1.0;
w2=w/2.0;
dh=w/400.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  w,0.0,dx};
Point(4)={0.0,  w,0.0,dx};

Point(5)={0.0,w2+dh,0.0,dx};
Point(6)={ w2,w2+dh,0.0,dx};
Point(7)={ w2,   w2,0.0,dx};
Point(8)={ w2,w2-dh,0.0,dx};
Point(9)={0.0,w2-dh,0.0,dx};//+
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
Line(6) = {8, 9};
//+
Line(7) = {9, 1};
//+
Circle(8) = {8, 7, 6};
//+
Line Loop(1) = {4, 5, -8, 6, 7, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {4, 7};
//+
Physical Line("bottom") = {1};
//+
Physical Line("right") = {2};
//+
Physical Line("top") = {3};
//+
Physical Surface("block") = {1};
