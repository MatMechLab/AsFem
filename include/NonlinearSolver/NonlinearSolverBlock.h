//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_NONLINEARSOLVERINFO_H
#define ASFEM_NONLINEARSOLVERINFO_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

#include "NonlinearSolverType.h"

using namespace std;

class NonlinearSolverBlock{
public:
    string _SolverName;
    NonlinearSolverType _SolverType;
    LineSearchType _LineSearchType;
    PetscInt _MaxIters=25,_LineSearchOrder;
    PetscReal _RAbsTol=1.0e-8,_RRelTol=1.0e-10,_STol=1.0e-16;

    void Reset(){
        // _SolverName="NewtonRaphson";
        // _SolverType=NonlinearSolverType::NewtonRaphson;
        // the following has better convergence behavior
        _SolverName="Newton line search";
        _SolverType=NonlinearSolverType::SNESNewtonLs;
        _LineSearchType=LineSearchType::LineSearchDefault;
        _LineSearchOrder=2;
        _MaxIters=25;
        _RAbsTol=1.0e-8;
        _RRelTol=1.0e-10;
        _STol=1.0e-16;// |dx|<|x|*stol
    }

    void PrintNonlinearSolverBlock()const{
        PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
        PetscPrintf(PETSC_COMM_WORLD,"*** Nonlinear solver block information:                               ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***   nonlinear solver name = [%38s]***\n",_SolverName.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"***   max iterations        = %4d, stol=%14.6e               ***\n",_MaxIters,_STol);
        PetscPrintf(PETSC_COMM_WORLD,"***   residual rel error    = %14.6e,abs error=%14.6e ***\n",_RRelTol,_RAbsTol);
        if(_LineSearchType==LineSearchType::LineSearchDefault){
            PetscPrintf(PETSC_COMM_WORLD,"***   line search type      = default       , order= %2d               ***\n",_LineSearchOrder);
        }
        else if(_LineSearchType==LineSearchType::LineSearchBasic){
            PetscPrintf(PETSC_COMM_WORLD,"***   line search type      = basic         , order= %2d               ***\n",_LineSearchOrder);
        }
        else if(_LineSearchType==LineSearchType::LineSearchBackTrace){
            PetscPrintf(PETSC_COMM_WORLD,"***   line search type      = backtrace     , order= %2d               ***\n",_LineSearchOrder);
        }
        else if(_LineSearchType==LineSearchType::LineSearchL2){
            PetscPrintf(PETSC_COMM_WORLD,"***   line search type      = l2 norm       , order= %2d               ***\n",_LineSearchOrder);
        }
        else if(_LineSearchType==LineSearchType::LineSearchCP){
            PetscPrintf(PETSC_COMM_WORLD,"***   line search type      = critical point, order= %2d               ***\n",_LineSearchOrder);
        }
    }
};

#endif //ASFEM_NONLINEARSOLVERINFO_H