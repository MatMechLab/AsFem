dx=0.5;

r1=60.0;
r2=5.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={ r1,0.0,0.0,dx};
Point(3)={0.0, r1,0.0,dx};


Point(4)={r1+r2,  0.0,0.0,dx};
Point(5)={  0.0,r1+r2,0.0,dx};//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 5};
//+
Circle(5) = {3, 1, 2};
//+
Circle(6) = {5, 1, 4};
//+
Line Loop(1) = {3, 5, -1};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {6, -2, -5, 4};
//+
Plane Surface(2) = {2};
//+
Physical Surface("inner") = {1};
//+
Physical Surface("outer") = {2};
