dx=0.1;

w=10.0;h=10.0;

Point(1)={  0.0,  0.0,0.0,dx};
Point(2)={w/2.0,  0.0,0.0,dx};
Point(3)={w/2.0,h/2.0,0.0,dx};
Point(4)={    w,h/2.0,0.0,dx};
Point(5)={    w,    h,0.0,dx};
Point(6)={  0.0,    h,0.0,dx};//+
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
Line(6) = {6, 1};
//+
Line Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {6};
//+
Physical Line("bottom1") = {1};
//+
Physical Line("bottom2") = {3};
//+
Physical Line("right1") = {2};
//+
Physical Line("right2") = {4};
//+
Physical Line("top") = {5};
//+
Physical Surface("block") = {1};
