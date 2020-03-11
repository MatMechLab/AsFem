//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "NonlinearSolver/NonlinearSolver.h"

NonlinearSolver::NonlinearSolver(){
    _Rnorm0=1.0;_Rnorm=1.0;
    _Enorm0=1.0;_Enorm=1.0;
    _RAbsTol=1.0e-8;_RRelTol=1.0e-10;
    _EAbsTol=1.0e-19;_ERelTol=1.0e-20;
    _MaxIters=100;_Iters=0;
    _STol=1.0e-20;
    _SolverType=NonlinearSolverType::NewtonRaphson;
}

//***************************************
// PetscErrorCode FormJacobian(SNES snes,Vec x,Mat A,Mat B,void *ctx){
//     AppCtx *user=(AppCtx*)ctx;
//     VecCopy(user->_solution._Uold,x);
//     user->_feSystem.FormFE(6,1.0,1.0,ctan,
//     user->_mesh,user->_dofHandler,user->_fe,user->_elmtSystem,
//     user->_mateSystem,x,user->)
// }
// PetscErrorCode FormResidual(SNES snes,Vec x,Vec RHS,void *ctx);