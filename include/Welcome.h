//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "petsc.h"


void Welcome(const PetscInt &year,const PetscInt &month,const PetscInt &day,const PetscReal &version){
    PetscInt Major,Minor,SubMinor;
    PetscGetVersionNumber(&Major,&Minor,&SubMinor,NULL);
    PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Welcome to use AsFem                                              ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** A Simple Finite Element Method program                            ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Version: %-10.2f  Release @ %4d-%02d-%02d                         ***\n",version,year,month,day);
    PetscPrintf(PETSC_COMM_WORLD,"*** PETSc version: %2d.%2d.%-2d                                           ***\n",Major,Minor,SubMinor);
    PetscPrintf(PETSC_COMM_WORLD,"*** License: GPL-3.0                                                  ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Author: Yang Bai                                                  ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Contact: walkandthinker@gmail.com                                 ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** QQ Group: 879908352                                               ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Website: https://github.com/yangbai90/AsFem                       ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Feel free to use and discuss  .:.                                 ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
}
