//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.26
//+++ Purpose: Define equation system in AsFem, here you can access
//+++          K matrix and Residual of our system equations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "EquationSystem/EquationSystem.h"

EquationSystem::EquationSystem(){
    _nDofs=0;
}
//**************************************************
void EquationSystem::InitEquationSystem(const int &ndofs,const int &maxrownnz){
    _nDofs=ndofs;

    VecCreate(PETSC_COMM_WORLD,&_RHS);
    VecSetSizes(_RHS,PETSC_DECIDE,_nDofs);
    VecSetUp(_RHS);
    VecSet(_RHS,0.0);

    //***************************************************************
    //*** here the maxrownnz should come from our dofhandler, where we create the dof map,
    //*** thereby, we can get the maximum non-zero entities of all the rows
    //***************************************************************
    MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,_nDofs,_nDofs,maxrownnz,NULL,maxrownnz,NULL,&_AMATRIX);

    //*************************************************************************************************************
    //*** here we allow PETSc to allocate or extend the width of each row in our matrix
    //*** after the first element-loop, once our matrix is initialized, we should disable the new allocation !!!
    //*************************************************************************************************************
    MatSetOption(_AMATRIX,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
}
//*********************************************************************************

void EquationSystem::ReleaseMem(){
    MatDestroy(&_AMATRIX);
    VecDestroy(&_RHS);
}