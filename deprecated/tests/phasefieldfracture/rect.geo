dx=0.04;

n1=51;
n2=81;

w=1.0;h=1.0;
w2=w/2.0;
h2=h/2.0;
dh=h/250.0;


Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  h,0.0,dx};
Point(4)={0.0,  h,0.0,dx};

Point(5)={w/2.0,h/2.0,0.0,dx};

Point(6)={0.0,h2+dh,0.0,dx};
Point(7)={0.0,h2-dh,0.0,dx};

Point(8)={w,h/2.0,0.0,dx};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 8};
//+
Line(3) = {8, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {5, 7};
//+
Line(8) = {7, 1};
//+
Line Loop(1) = {5, 6, 7, 8, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {5, 8};
//+
Physical Line("right") = {3, 2};
//+
Physical Line("bottom") = {1};
//+
Physical Line("top") = {4};
//+
Physical Surface("block") = {1};
//+
Transfinite Line {5, 6, 7, 8, 3, 2} = n1 Using Progression 1;
//+
Transfinite Line {4, 1} = n2 Using Progression 1;
