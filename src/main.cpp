//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include <iostream>
#include "petsc.h"

#include "Welcome.h"

#include "FEProblem/FEProblem.h"

int main(int args,char *argv[]){
    PetscErrorCode ierr;
    ierr=PetscInitialize(&args,&argv,NULL,NULL);if (ierr) return ierr;

    const PetscInt Year=2020;
    const PetscInt Month=2;
    const PetscInt Day=14;
    const PetscReal Version=0.5;

    Welcome(Year,Month,Day,Version);
    
    FEProblem feProblem(args,argv);

    feProblem.RunFEProblem();


    ierr=PetscFinalize();CHKERRQ(ierr);
    return ierr;
}