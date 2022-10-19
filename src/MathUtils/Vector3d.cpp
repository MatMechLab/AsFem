//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.13
//+++ Purpose: defines the vector with only 3-components, this will
//+++          be frequently used in shape function calculation.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/Vector3d.h"

Vector3d::Vector3d(){
    m_vals[0]=0.0;m_vals[1]=0.0;m_vals[2]=0.0;
}
Vector3d::Vector3d(const double &val){
    m_vals[0]=val;m_vals[1]=val;m_vals[2]=val;
}
Vector3d::Vector3d(const Vector3d &a){
    m_vals[0]=a.m_vals[0];m_vals[1]=a.m_vals[1];m_vals[2]=a.m_vals[2];
}

Vector3d operator*(const double &val,const Vector3d &a){
    Vector3d temp(0.0);
    temp(1)=val*a(1);temp(2)=val*a(2);temp(3)=val*a(3);
    return temp;
}