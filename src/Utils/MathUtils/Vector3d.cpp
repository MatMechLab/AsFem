//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.10.17
//+++ Purpose: Implement the vector (3 components) operator and 
//+++          calculation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/Vector3d.h"

Vector3d::Vector3d(const double &val){
    _vals[0]=val;_vals[1]=val;_vals[2]=val;
}

Vector3d operator*(const double &val,const Vector3d &a){
    Vector3d temp(0.0);
    temp(1)=a(1)*val;
    temp(2)=a(2)*val;
    temp(3)=a(3)*val;
    return temp;
}

ostream& operator<<(ostream& os,const Vector3d &vec){
    os<<vec(1)<<" ,"<<vec(2)<<" ,"<<vec(3)<<endl;
    return os;
}