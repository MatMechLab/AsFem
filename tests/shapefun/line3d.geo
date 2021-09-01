SetFactory("OpenCASCADE");

dx=0.25;
x=1.0;y=2.0;z=3.0;

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
