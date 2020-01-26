dx=0.05;

R=2.5;w=2.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  R,0.0,0.0,dx};
Point(3)={R+w,0.0,0.0,dx};
Point(4)={0.0,R+w,0.0,dx};
Point(5)={0.0,  R,0.0,dx};

//+
Circle(1) = {5, 1, 2};
//+
Circle(2) = {4, 1, 3};
//+
Line(3) = {4, 5};
//+
Line(4) = {3, 2};
//+
Line Loop(1) = {1, -4, -2, 3};
//+
Plane Surface(1) = {1};

//+
Physical Line("left") = {3};
//+
Physical Line("bottom") = {4};
//+
Physical Line("surfacei") = {1};
//+
Physical Line("surfaceo") = {2};

//+
Physical Surface("block") = {1};
