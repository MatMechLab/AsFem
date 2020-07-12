//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: define the nonlinear solver class in AsFem
//+++          this class mainly call the SNES subroutines of PETSc
//+++          to solve the nonlinear equations
//+++          once PETSc update the API, we should also update the
//+++          related code !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

#include "NonlinearSolver/NonlinearSolverBlock.h"

class NonlinearSolver{
public:
    NonlinearSolver();
    void Init(NonlinearSolverBlock nonlinearsolverblock);

private:
    //*********************************************
    //*** For nonlinear solver information
    //*********************************************
    double _Rnorm0,_Rnorm;
    double _Enorm0,_Enorm;
    double _dUnorm,_dUnorm0;
    double _RAbsTol,_RRelTol;
    double _EAbsTol,_ERelTol;
    double _STol;
    int _MaxIters,_Iters;
    bool _IsConvergent;
    NonlinearSolverType _SolverType;
    string _SolverName;
    string _PCTypeName;
};