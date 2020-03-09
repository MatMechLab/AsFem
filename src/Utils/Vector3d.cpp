//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/Vector3d.h"

inline Vector3d operator*(const double &val,const Vector3d &a){
    Vector3d temp(0.0);
    temp._vals[0]=a._vals[0]*val;
    temp._vals[1]=a._vals[1]*val;
    temp._vals[2]=a._vals[2]*val;
    return temp;
}