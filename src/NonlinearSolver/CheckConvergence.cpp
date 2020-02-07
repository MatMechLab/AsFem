//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2019
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "NonlinearSolver/NonlinearSolver.h"

bool NonlinearSolver::CheckConvergence(){
    if(_Rnorm<_RAbsTol||_Rnorm<_Rnorm0*_RRelTol){
        return true;
    }
    else if(_Enorm<_EAbsTol||_Enorm<_Enorm0*_ERelTol){
        return true;
    }
    return false;
}