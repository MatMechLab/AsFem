dx=0.05;

w=5.0;r=1.0;
n1=81;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={w-r,0.0,0.0,dx};
Point(3)={  w,0.0,0.0,dx};
Point(4)={  w,  r,0.0,dx};
Point(5)={  w,  w,0.0,dx};
Point(6)={0.0,  w,0.0,dx};//+
Line(1) = {1, 2};
//+
Line(2) = {4, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 1};
//+
Circle(5) = {2, 3, 4};
//+
Line Loop(1) = {4, 1, 5, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {4};
//+
Physical Line("bottom") = {1};
//+
Physical Line("right") = {2};
//+
Physical Line("top") = {3};
//+
Physical Line("surface") = {5};
//+
Physical Surface("block") = {1};
//+
Transfinite Line {5} = n1 Using Progression 1;
