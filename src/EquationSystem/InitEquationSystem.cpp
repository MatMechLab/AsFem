//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "EquationSystem/EquationSystem.h"

void EquationSystem::InitEquationSystem(const PetscInt &ndofs,const PetscInt &maxrownnz){
    _nDofs=ndofs;

    VecCreate(PETSC_COMM_WORLD,&_RHS);
    VecSetSizes(_RHS,PETSC_DECIDE,_nDofs);
    VecSetUp(_RHS);
    VecSet(_RHS,0.0);

    MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,_nDofs,_nDofs,maxrownnz,NULL,maxrownnz,NULL,&_AMATRIX);
    MatSetOption(_AMATRIX,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
    //MatSetUp(_AMATRIX);
}