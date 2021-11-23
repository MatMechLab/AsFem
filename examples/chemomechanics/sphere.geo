SetFactory("OpenCASCADE");

dx=0.025;

r=1.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  r,0.0,0.0,dx};
Point(3)={0.0,  r,0.0,dx};//+
Line(1) = {1, 2};
//+
Line(2) = {1, 3};
//+
Circle(3) = {3, 1, 2};
//+
Curve Loop(1) = {2, 3, -1};
//+
Plane Surface(1) = {1};
//+
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{1}; 
}
//+
Physical Surface("ux", 8) = {4};
//+
Physical Surface("uy", 9) = {3};
//+
Physical Surface("uz", 10) = {1};
//+
Physical Surface("surface", 11) = {2};
//+
Physical Volume("block", 12) = {1};
