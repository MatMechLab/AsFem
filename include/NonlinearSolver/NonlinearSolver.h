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
//+++ Date   : 2022.08.20
//+++ Purpose: the nonlinear solver class in AsFem, users should
//+++          offer the details by inheriting the nlsolver base class.
//+++          this class should inherit all the built-in and 
//+++          user-defined nonlinear solvers for R(x)->0 problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "NewtonRaphsonSolver.h"
#include "LinearSolver/LinearSolver.h"
#include "NonlinearSolver/SNESSolver.h"


/**
 * This class implement and manage all the nonlinear solvers in AsFem.
 * The R(x)->0 problem will be solved within this class
 */
class NonlinearSolver:public SNESSolver,
                      public NewtonRaphsonSolver{
public:
    /**
     * constructor
     */
    NonlinearSolver();

    /**
     * init the nonlinear solver
     * @param lsolver the linear solver class
     */
    void init(LinearSolver &lsolver);

    /**
     * solve the nonlinear equation, if success then return true
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dof class
     * @param t_FE the fe class
     * @param t_ElmtSystem the element system class
     * @param t_MateSystem the material system class
     * @param t_FESystem the fe system class
     * @param t_BCSystem the boundary condition system
     * @param t_SolnSystem the solution system class
     * @param t_EqSystem the equation system class
     * @param t_LinearSolver the linear solver system
     * @param t_FECtrlInfo the fe control info
     */
    virtual bool solve(FECell &t_FECell,
                       DofHandler &t_DofHandler,
                       FE &t_FE,
                       ElmtSystem &t_ElmtSystem,
                       MateSystem &t_MateSystem,
                       FESystem &t_FESystem,
                       BCSystem &t_BCSystem,
                       SolutionSystem &t_SolnSystem,
                       EquationSystem &t_EqSystem,
                       LinearSolver &t_LinearSolver,
                       FEControlInfo &t_FECtrlInfo) final;

    /**
     * get the nonlinear iterations
     * @return return the nonlinear iterations
     */
    int getIterationNum()const;

    void releaseMemory();

    void printSolverInfo()const;

public:
    NonlinearSolverBlock m_NlSolverBlock;/**< the nonlinear solver block defined in json file */

};