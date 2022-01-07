//+
SetFactory("OpenCASCADE");

dx=0.5;
w=1.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={  w,0.0,0.0,dx};

Line(1)={1,2};//+
Physical Point("left") = {1};
//+
Physical Curve("block") = {1};
