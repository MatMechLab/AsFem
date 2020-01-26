#ifndef ASFEM_POINT_H
#define ASFEM_POINT_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

class Point
{
public:
    Point() {_coords[0]=0.0;_coords[1]=0.0;_coords[2]=0.0;}
    Point(const double &x,double &y,double &z) {_coords[0]=x;_coords[1]=y;_coords[2]=z;}

    inline double  operator()(const int &i) const {return _coords[i-1];}// must start from 1
    inline double& operator()(const int &i) {return _coords[i-1];}

private:
    double _coords[3];
};

#endif // ASFEM_POINT_H