//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "NonlinearSolver/NonlinearSolver.h"

void NonlinearSolver::PrintIterationInfo()const{
    PetscPrintf(PETSC_COMM_WORLD,"***  Iters=%3d,|R0|=%12.5e,|E0|=%12.5e,|dU0|=%12.5e ***\n",_Iters,_Rnorm0,_Enorm0,_dUnorm0);
    PetscPrintf(PETSC_COMM_WORLD,"***            |R |=%12.5e,|E |=%12.5e,|dU |=%12.5e ***\n",_Rnorm,_Enorm,_dUnorm);
}
//*****************************
void NonlinearSolver::PrintIterationDetailsInfo()const{
    PetscPrintf(PETSC_COMM_WORLD,"***  Iters=%3d,|R|=%13.5e,|E|=%13.5e,|dU|=%13.5e ***\n",_Iters,_Rnorm,_Enorm,_dUnorm);
}