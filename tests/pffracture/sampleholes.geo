SetFactory("OpenCASCADE");

dx=0.01;
w=1.0;

n=81;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  w,0.0,dx};
Point(4)={0.0,  w,0.0,dx};//+

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};

Circle(5) = {0.25, 0.15, 0, 0.08, 0, 2*Pi};
//+
Circle(6) = {0.5, 0.15, 0, 0.08, 0, 2*Pi};
//+
Circle(7) = {0.75, 0.15, 0, 0.08, 0, 2*Pi};

//+
Circle(8) = {0.25, 0.4, 0, 0.08, 0, 2*Pi};
//+
Circle(9) = {0.5, 0.4, 0, 0.08, 0, 2*Pi};
//+
Circle(10) = {0.75, 0.4, 0, 0.08, 0, 2*Pi};

//+
Circle(11) = {0.75, 0.625, 0, 0.08, 0, 2*Pi};
//+
Circle(12) = {0.5, 0.625, 0, 0.08, 0, 2*Pi};
//+
Circle(13) = {0.25, 0.625, 0, 0.08, 0, 2*Pi};
//+

Circle(14) = {0.25, 0.85, 0, 0.08, 0, 2*Pi};
//+
Circle(15) = {0.5, 0.85, 0, 0.08, 0, 2*Pi};
//+
Circle(16) = {0.75, 0.85, 0, 0.08, 0, 2*Pi};
//+
Line Loop(2) = {14};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {15};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {16};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {13};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {12};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {11};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {8};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {9};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {10};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {5};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {6};
//+
Plane Surface(12) = {12};
//+
Line Loop(13) = {7};
//+
Plane Surface(13) = {13};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }
//+
Physical Line("left") = {17};
//+
Physical Line("bottom") = {19};
//+
Physical Line("right") = {20};
//+
Physical Line("top") = {18};
//+
Physical Surface("block") = {1};
//+
Transfinite Line {14, 15, 16, 13, 12, 11, 8, 9, 10, 5, 6, 7} = n Using Progression 1;
