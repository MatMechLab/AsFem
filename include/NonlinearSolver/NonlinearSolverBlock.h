//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
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


using std::vector;
using std::string;

/**
 * this class stores the basic information of [nonlinearsolver] block
 */
class NonlinearSolverBlock{
public:
    NonlinearSolverBlock(){
        m_NlSolverTypeName="newton with line search";
        m_NlSolverType=NonlinearSolverType::NEWTONLS;
        m_MaxIters=25;
        m_AbsTolR=7.5e-7;
        m_RelTolR=5.0e-10;
        m_STol=0.0; // |dx|<|x|*stol
        m_PCTypeName="lu";
        m_LinearSolverName="default(gmres)";
        m_CheckJacobian=false;
    }

    string              m_NlSolverTypeName;/**< the string name of nonlinear solver */
    NonlinearSolverType m_NlSolverType; /**< the type of nonlinear solver */
    int                 m_MaxIters;/**< maximum iterations */
    double m_AbsTolR;/**< the absolute tolerance for residual */
    double m_RelTolR;/**< the relative tolerance for residual */
    double m_STol;/**< the tolerance for line search */
    string m_LinearSolverName;/**< the linear solver name */

    string m_PCTypeName;/**< the string name of preconditioner */
    bool m_CheckJacobian=false;/**< if this is true, then SNES will compare your jacobian with the finite difference one */

    /**
     * initialize the nlsolver block
     */
    void init(){
        m_NlSolverTypeName="newton with line search";
        m_NlSolverType=NonlinearSolverType::NEWTONLS;
        m_MaxIters=25;
        m_AbsTolR=7.5e-7;
        m_RelTolR=5.0e-10;
        m_STol=0.0; // |dx|<|x|*stol
        m_PCTypeName="lu";
        m_LinearSolverName="default(gmres)";
        m_CheckJacobian=false;
    }
};