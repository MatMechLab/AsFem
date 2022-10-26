//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
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


using std::vector;
using std::string;

/**
 * this class stores the basic information of [nonlinearsolver] block
 */
class NonlinearSolverBlock{
public:
    NonlinearSolverBlock(){
        m_nlsolvertypename="newton with line search";
        m_nlsolvertype=NonlinearSolverType::NEWTONLS;
        m_maxiters=25;
        m_abstol_r=7.5e-7;
        m_reltol_r=5.0e-10;
        m_s_tol=0.0; // |dx|<|x|*stol
        m_pctypename="lu";
        m_linearsolvername="default(gmres)";
        m_checkjacobian=false;
    }

    string              m_nlsolvertypename;/**< the string name of nonlinear solver */
    NonlinearSolverType m_nlsolvertype; /**< the type of nonlinear solver */
    int                 m_maxiters;/**< maximum iterations */
    double m_abstol_r;/**< the absolute tolerance for residual */
    double m_reltol_r;/**< the relative tolerance for residual */
    double m_s_tol;/**< the tolerance for line search */
    string m_linearsolvername;/**< the linear solver name */

    string m_pctypename;/**< the string name of preconditioner */
    bool m_checkjacobian=false;/**< if this is true, then SNES will compare your jacobian with the finite difference one */

    /**
     * initialize the nlsolver block
     */
    void init(){
        m_nlsolvertypename="newton with line search";
        m_nlsolvertype=NonlinearSolverType::NEWTONLS;
        m_maxiters=25;
        m_abstol_r=7.5e-7;
        m_reltol_r=5.0e-10;
        m_s_tol=0.0; // |dx|<|x|*stol
        m_pctypename="lu";
        m_linearsolvername="default(gmres)";
        m_checkjacobian=false;
    }
};