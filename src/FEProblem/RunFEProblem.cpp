//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FEProblem/FEProblem.h"

void FEProblem::RunFEProblem(){
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    if(_rank==0){
        _JobStartTime=chrono::system_clock::now();
        auto in_time_t=chrono::system_clock::to_time_t(_JobStartTime);
        stringstream ss;
        ss<<put_time(localtime(&in_time_t),"%Y-%m-%d %X");
        _JobStartTimeContext=ss.str();
    }

    // PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** AsFem job starts at                     %25s ***\n",_JobStartTimeContext.c_str());
    // PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
    PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Start to read the input file ...                                  ***\n");
        
    if(!RunInputSystem()){
        PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: read input file failed due to some errors !!!              ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
        Msg_AsFem_Exit();
    }
    PetscPrintf(PETSC_COMM_WORLD,"*** Read input file finished !                                        ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
    
    InitFEProblem();
    RunAnalysis();


    //*************************************
    //*** job finished!
    //*************************************
    if(_rank==0){
        _JobEndTime=chrono::system_clock::now();
        auto in_time_t=chrono::system_clock::to_time_t(_JobEndTime);
        stringstream ss;
        ss<<put_time(localtime(&in_time_t),"%Y-%m-%d %X");
        _JobEndTimeContext=ss.str();
    }
    PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** AsFem job finished at                   %25s ***\n",_JobEndTimeContext.c_str());
    PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
}