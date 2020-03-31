//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FEProblem/FEProblem.h"

void FEProblem::RunAnalysis(){
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    if(_rank==0){
        _TimerStartOfJob=chrono::high_resolution_clock::now();
    }

    PetscPrintf(PETSC_COMM_WORLD,"*** Start to run FEM analysis ...                                     ***\n");
    if(_jobBlock._JobType==JobType::StaticJob){
        RunStaticAnalysis();
    }
    else if(_jobBlock._JobType==JobType::TransientJob){
        RunTransientAnalysis();
    }

    if(_rank==0){
        _TimerEndOfJob=chrono::high_resolution_clock::now();
        _DurationOfJob=Duration(_TimerStartOfJob,_TimerEndOfJob);
    }
    PetscPrintf(PETSC_COMM_WORLD,"*** Analysis finished                         !!!===>[%13.5e s]***\n",_DurationOfJob);
    
}