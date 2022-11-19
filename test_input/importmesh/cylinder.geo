//+
SetFactory("OpenCASCADE");

dx=1.0;
h=1.0;
r=0.5;

Point(1)={0.0,0.0,0.0,dx};//+
Circle(1) = {0, 0, 0, r, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 1} {
  Surface{1}; 
}
//+
Physical Surface("bottom", 4) = {1};
//+
Physical Surface("top", 5) = {3};
//+
Physical Surface("surface", 6) = {2};
//+
Physical Volume("block", 7) = {1};
