//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "NonlinearSolver/NonlinearSolver.h"

bool NonlinearSolver::LinearSolve(Mat &A,Vec &x,Vec &F){
    KSPSetOperators(_ksp,A,A);
    KSPSolve(_ksp,F,x);
    KSPGetConvergedReason(_ksp,&_reason);
    if(_reason==KSP_CONVERGED_RTOL){
        return true;
    }
    else if(_reason==KSP_CONVERGED_ATOL){
        return true;
    }
    else if(_reason==KSP_CONVERGED_ITS){
        return true;
    }
    else{
        return false;
    }
}