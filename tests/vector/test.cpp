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
//+++ Date   : 2022.01.15
//+++ Purpose: test cpp for the functionality of AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "petsc.h"

#include "Utils/Vector3d.h"

int main(int args,char *argv[]){
    PetscErrorCode ierr;
    ierr=PetscInitialize(&args,&argv,NULL,NULL);if (ierr) return ierr;

    Vector3d myvec;
    
    myvec(1)=numeric_limits<double>::epsilon();
    myvec(2)=numeric_limits<double>::epsilon();
    myvec(3)=numeric_limits<double>::epsilon();

    cout<<"norm is:"<<myvec.norm()<<", normsq is:"<<myvec.normsq()<<endl;

    myvec(1)=numeric_limits<long double>::epsilon();
    myvec(2)=numeric_limits<long double>::epsilon();
    myvec(3)=numeric_limits<long double>::epsilon();

    cout<<"norm is:"<<myvec.norm()<<", normsq is:"<<myvec.normsq()<<endl;

    ierr=PetscFinalize();CHKERRQ(ierr);
    return ierr;
}
