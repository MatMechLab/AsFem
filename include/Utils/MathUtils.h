//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_MATHUTILS_H
#define ASFEM_MATHUTILS_H

#include <iostream>

// #include "petsc.h"

// #include "Vector3d.h"

using namespace std;

// inline Vector3d operator*(const double &a,const Vector3d &b){
//     Vector3d temp;
//     temp(1)=a*b(1);temp(2)=a*b(2);temp(3)=a*b(3);
//     return temp;
// }

// inline double operator*(const Vector3d &a,const Vector3d &b){
//     return a(1)*b(1)+a(2)*b(2)+a(3)*b(3);
// }

inline double sign(const double &x){
    return x >= 0.0 ? 1.0 : -1.0;
}

// double AtDotB(const Vec &a,const Vec &b){
//     PetscScalar res;
//     VecTDot(a,b,&res);
//     return res;
// }

// Vec operator+(const Vec &a,const Vec &b){
//     Vec vnew;
//     VecDuplicate(a,&vnew);
//     VecWAXPY(vnew,1.0,a,b);
//     return vnew;
// }
// inline Vec& operator-(const Vec &a,const Vec &b){
//     Vec vnew;
//     VecDuplicate(a,&vnew);
//     VecCopy(a,vnew);//vnew=a
//     VecAXPY(vnew,-1.0,b);//vnew=vnew-b
//     return vnew;
// }
// inline Vec& operator=(const Vec &a){
//     Vec vnew;
//     VecDuplicate(a,&vnew);
//     VecCopy(a,vnew);
//     return vnew;
// }


#endif // ASFEM_MATHUTILS_H