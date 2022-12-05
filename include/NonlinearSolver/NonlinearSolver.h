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
//+++ Date   : 2022.08.20
//+++ Purpose: the nonlinear solver class in AsFem, users should
//+++          offer the details by inheriting the nlsolver base class.
//+++          this class should inherit all the built-in and 
//+++          user-defined nonlinear solvers for R(x)->0 problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include "NonlinearSolver/SNESSolver.h"


/**
 * This class implement and manage all the nonlinear solvers in AsFem.
 * The R(x)->0 problem will be solved within this class
 */
class NonlinearSolver:public SNESSolver{
public:
    /**
     * constructor
     */
    NonlinearSolver();

    /**
     * init the nonlinear solver
     */
    void init();

public:
    NonlinearSolverBlock m_nlsolverblock;/**< the nonlinear solver block defined in json file */

};