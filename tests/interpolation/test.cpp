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

#include "Utils/MathFuns.h"

int main(int args,char *argv[]){
    PetscErrorCode ierr;
    ierr=PetscInitialize(&args,&argv,NULL,NULL);if (ierr) return ierr;

    vector<double> xs,ys,xi;
    
    xs.push_back(0.0);ys.push_back(0.0);
    xs.push_back(1.0);ys.push_back(1.0);
    xs.push_back(2.0);ys.push_back(-1.0);

    xi.push_back(0.5);
    xi.push_back(1.5);
    xi.push_back(1.8);
    xi.push_back(2.5);
    xi.push_back(3.5);

    for(const auto &it:xi){
        cout<<"xi="<<it<<", interpolated value="<<PicewiseLinearInterpolation(xs,ys,it)<<endl;
    }

    ierr=PetscFinalize();CHKERRQ(ierr);
    return ierr;
}
