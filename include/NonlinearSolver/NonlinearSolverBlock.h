//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: define the [nonlinearsolver] block for our input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include "NonlinearSolver/NonlinearSolverType.h"


using namespace std;

/**
 * this class stores the basic information of [nonlinearsolver] block
 */
class NonlinearSolverBlock{
public:
    NonlinearSolverBlock(){
        _SolverTypeName="newton with line search";
        _SolverType=NonlinearSolverType::NEWTONLS;
        _MaxIters=25;
        _RAbsTol=4.5e-8;
        _RRelTol=1.0e-9;
        _STol=0.0; // |dx|<|x|*stol
        _PCTypeName="lu";
        _LinearSolverName="gmres";
        _CheckJacobian=false;
    }

    string              _SolverTypeName;
    NonlinearSolverType _SolverType;
    int                 _MaxIters;
    double _RAbsTol,_RRelTol,_STol;
    string _LinearSolverName;// for linear solver, i.e., ksp, mumps, superlu_dist

    string _PCTypeName;
    bool _CheckJacobian=false;/**< if this is true, then SNES will compare your jacobian with the finite difference one */

    void Init(){
        _SolverTypeName="newton with line search";
        _SolverType=NonlinearSolverType::NEWTONLS;
        _MaxIters=25;
        _RAbsTol=4.5e-8;
        _RRelTol=1.0e-9;
        _STol=0.0; // |dx|<|x|*stol
        _PCTypeName="lu";
        _LinearSolverName="gmres";
        _CheckJacobian=false;
    }
};
