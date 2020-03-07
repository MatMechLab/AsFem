//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2019
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_JOBBLOCK_H
#define ASFEM_JOBBLOCK_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

#include "JobType.h"

using namespace std;

class JobBlock{
public:
    string _JobTypeName;
    JobType _JobType;
    bool _IsDebug,_IsDepDebug;
    bool _IsProjection;
    PetscInt _Interval;

    void Reset(){
        _JobTypeName.clear();
        _JobType=JobType::NullJob;
        _IsDebug=true;
        _IsDepDebug=false;
        _IsProjection=false;
        _Interval=0;
    }

    void PrintJobBlockInfo() const{
        PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
        PetscPrintf(PETSC_COMM_WORLD,"*** Job block information:                                            ***\n");
        
        if(_JobType==JobType::StaticJob){
            PetscPrintf(PETSC_COMM_WORLD,"***   job type        = static                                        ***\n");
        }
        else if(_JobType==JobType::TransientJob){
            PetscPrintf(PETSC_COMM_WORLD,"***   job type        = transient                                     ***\n");
        }

        PetscPrintf(PETSC_COMM_WORLD,"***   output interval = %2d                                            ***\n",_Interval);
        if(_IsDebug){
            if(_IsDepDebug){
                PetscPrintf(PETSC_COMM_WORLD,"***   debug           = deep                                          ***\n");
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"***   debug           = true                                          ***\n");
            }
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"***   debug           = false                                          ***\n");   
        }
        if(_IsProjection){
            PetscPrintf(PETSC_COMM_WORLD,"***   projection      = true                                          ***\n");
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"***   projection      = false                                         ***\n");
        }
        PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");

    }
};


#endif // ASFEM_JOBBLOCK_H