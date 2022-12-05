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
//+++ Date   : 2020.12.28
//+++ Purpose: the main program of the whole AsFem framework
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "petsc.h"

#include "Welcome.h"
#include "FEProblem/FEProblem.h"
#include "MathUtils/SparseMatrix.h"

int main(int args,char *argv[]){
    PetscErrorCode ierr;
    ierr=PetscInitialize(&args,&argv,NULL,NULL);if (ierr) return ierr;

    const PetscInt Year=2022;
    const PetscInt Month=10;
    const PetscInt Day=19;
    const PetscReal Version=0.8;

    welcome(Year,Month,Day,Version);

    FEProblem feProblem;
    feProblem.initFEProblem(args,argv);
    feProblem.run();
    feProblem.finalize();
   

    ierr=PetscFinalize();CHKERRQ(ierr);
    return ierr;
}
