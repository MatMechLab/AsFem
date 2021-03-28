//+
SetFactory("OpenCASCADE");

dx=0.25;
dh=5.0e-3;
dz=1.2;

n1=101;// two long notched edges
n2=31;// two short notched edges

Point(1)={0.0,0.0,0.0,dx};
Point(2)={5.0,0.0,0.0,dx};
Point(3)={5.0,5.0,0.0,dx};
Point(4)={5.0,10.,0.0,dx};
Point(5)={0.0,10.,0.0,dx};

Point(6)={5.0-2.0,5.0,0.0,dx};
Point(7)={5.0,5.0-dh,0.0,dx};
Point(8)={5.0,5.0+dh,0.0,dx};
Point(9)={0.0,5.0,0.0,dx};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 7};
//+
Line(3) = {7, 6};
//+
Line(4) = {6, 8};
//+
Line(5) = {8, 4};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 9};
//+
Line(8) = {9, 6};
//+
Line(9) = {9, 1};
//+
Curve Loop(1) = {7, 8, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 1, 2, 3, -8};
//+
Plane Surface(2) = {2};
//+
Extrude {0, 0, dz} {
  Surface{1}; Surface{2}; 
}
//+
Physical Surface("bottom") = {10};
//+
Physical Surface("top") = {7};
//+
Physical Volume("block") = {1, 2};
//+
Transfinite Curve {14, 8} = n1 Using Progression 1;
//+
Transfinite Curve {11, 13} = n2 Using Progression 1;
//+
Transfinite Surface {4};
