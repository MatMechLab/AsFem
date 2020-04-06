//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_TIMESTEPPINGBLOCK_H
#define ASFEM_TIMESTEPPINGBLOCK_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

#include "TimeSteppingType.h"

using namespace std;


class TimeSteppingBlock{
public:
    string _TimeSteppingMethodName;
    TimeSteppingType _TimeSteppingMethod=TimeSteppingType::BackWardEuler;
    PetscReal _FinalTime=1.0e-5;
    PetscReal _dtmax=0.1,_dtmin=1.0e-12,_dt0=1.0e-5;
    PetscInt _interval=1;
    PetscInt _nOpts=4;
    PetscReal _CutFactor=0.85,_GrowthFactor=1.1;
    bool _IsAdaptive=false;

    void Reset(){
        _TimeSteppingMethodName="BackwardEuler";
        _TimeSteppingMethod=TimeSteppingType::BackWardEuler;
        _FinalTime=1.0e-5;
        _dtmax=0.1;
        _dtmin=1.0e-12;
        _dt0=1.0e-5;
        _CutFactor=0.85;
        _GrowthFactor=1.1;
        _interval=1;
        _nOpts=4;
        _IsAdaptive=false;
    }

    void PrintTimeSteppingBlock()const{
        PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
        PetscPrintf(PETSC_COMM_WORLD,"*** Time stepping block information:                                  ***\n");
        if(_TimeSteppingMethod==TimeSteppingType::BackWardEuler){
            PetscPrintf(PETSC_COMM_WORLD,"***   time stepping method = Backward Euler                           ***\n");
        }
        else if(_TimeSteppingMethod==TimeSteppingType::CrankNicolson){
            PetscPrintf(PETSC_COMM_WORLD,"***   time stepping method = Crank Nicolson                           ***\n");
        }
        else if(_TimeSteppingMethod==TimeSteppingType::ThetaMethod){
            PetscPrintf(PETSC_COMM_WORLD,"***   time stepping method = Theta Method                             ***\n");
        }
        else if(_TimeSteppingMethod==TimeSteppingType::AlphaMethod){
            PetscPrintf(PETSC_COMM_WORLD,"***   time stepping method = Alpha Method                             ***\n");
        }
        PetscPrintf(PETSC_COMM_WORLD,"***   final time= %13.5e,    dt= %13.5e                 ***\n",_FinalTime,_dt0);
        PetscPrintf(PETSC_COMM_WORLD,"***   dtmax     = %13.5e, dtmin= %13.5e                 ***\n",_dtmax,_dtmin);
        
        if(_IsAdaptive){
            PetscPrintf(PETSC_COMM_WORLD,"***   output interval =%3d, adaptive=true, optimal iters=%3d          ***\n",_interval,_nOpts);
            PetscPrintf(PETSC_COMM_WORLD,"***   cut factor= %13.5e, growth factor=%13.5e          ***\n",_CutFactor,_GrowthFactor);
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"***   output interval =%3d, adaptive=false,optimal iters=%3d          ***\n",_interval,_nOpts);
        }
        // PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n"); 
    }
};

#endif // ASFEM_TIMESTEPPINGBLOCK_H