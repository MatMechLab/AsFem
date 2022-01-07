//+
SetFactory("OpenCASCADE");

dx=0.005;

r=0.085;

dw=0.025;

w=0.7;

n1=66;
n2=56;

//Mesh.MinimumCircleNodes=121;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};
Point(3)={  w,  w,0.0,dx};
Point(4)={0.0,  w,0.0,dx};//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Ellipse(5) = {0.175, 0.425, 0, 0.225, 0.175, 0, 2*Pi};
//+
Ellipse(6) = {0.6, 0.2, 0, 0.15, 0.1, 0, 2*Pi};
//+
Ellipse(7) = {0.825, 0.52, 0, 0.325, 0.1, 0, 2*Pi};
//+
Curve Loop(2) = {5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {6};
//+
Plane Surface(4) = {4};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Delete; }
//+
Extrude {0, 0, dw} {
  Surface{1}; 
}
//+
Physical Surface("left", 34) = {12, 9};
//+
Physical Surface("right", 35) = {3, 5, 7};
//+
Physical Surface("bottom", 36) = {8};
//+
Physical Surface("top", 37) = {2};
//+
Physical Surface("back", 38) = {1};
//+
Physical Surface("front", 39) = {13};
//+
Physical Volume("block", 40) = {1};
//+
Transfinite Curve {31, 27, 29, 15, 17, 19, 21} = 11 Using Progression 1;
//+
Transfinite Curve {33, 11, 20, 4} = 66 Using Progression 1;
//+
Transfinite Curve {8, 28} = 81 Using Progression 1;
//+
Transfinite Curve {32, 10, 30, 9} = 121 Using Progression 1;
//+
Transfinite Curve {18, 3} = 121 Using Progression 1;
//+
Transfinite Curve {22, 5} = 141 Using Progression 1;
//+
Transfinite Curve {16, 2} = 21 Using Progression 1;
//+
Transfinite Curve {24, 6} = 31 Using Progression 1;
