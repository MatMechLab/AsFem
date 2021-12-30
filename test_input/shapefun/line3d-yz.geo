SetFactory("OpenCASCADE");

dx=0.2;
x=0.0;y=1.0;z=2.0;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={x,y,z,dx};

//+
Line(1) = {1, 2};
//+
Physical Point("left", 2) = {1};
//+
Physical Point("right", 3) = {2};
//+
Physical Curve("block", 4) = {1};
