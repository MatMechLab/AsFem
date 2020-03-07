dx=0.025;

r1=4.0;
dr=1.5;


Point(1)={  0.0,  0.0,0.0,dx};
Point(2)={   r1,  0.0,0.0,dx};
Point(3)={r1+dr,  0.0,0.0,dx};
Point(4)={  0.0,   r1,0.0,dx};
Point(5)={  0.0,r1+dr,0.0,dx};
//+
Circle(1) = {4, 1, 2};
//+
Circle(2) = {5, 1, 3};
//+
Line(3) = {5, 4};
//+
Line(4) = {2, 3};
//+
Line Loop(1) = {3, 1, 4, -2};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {3};
//+
Physical Line("bottom") = {4};
//+
Physical Line("surface") = {1};
//+
Physical Line("surface2") = {2};
//+
Physical Surface("block") = {1};
