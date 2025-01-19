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
//+++ Date   : 2025.01.18
//+++ Purpose: the newton-raphson solver from AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "LinearSolver/LinearSolver.h"
#include "NonlinearSolver/NonlinearSolverBase.h"
#include "Utils/Timer.h"


class NewtonRaphsonSolver : public NonlinearSolverBase {
public:
    NewtonRaphsonSolver();

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
                       FEControlInfo &t_FECtrlInfo) override;

    inline int getNRIterationNum()const {
        return m_Iterations;
    }

    /**
     * set the maximum iterations of NR
     * @param t_MaxIterationNum the maximum nonolinear iterations
     */
    void setMaxIterationNum(int t_MaxIterationNum) {
        m_MaxIterations = t_MaxIterationNum;
    }

    /**
     * set the tolerance of the absolute error of residual
     * @param t_Tolerance the tolerance of residual
     */
    void setRAbsTolerance(double t_Tolerance) {
        m_RAbsTol=t_Tolerance;
    }

    /**
     * set the relative tolerance of the residual
     * @param t_Tolerance the tolerance
     */
    void setRRelTolerance(double t_Tolerance) {
        m_RRelTol=t_Tolerance;
    }

    void printSolverInfo()const;


private:
    int m_MaxIterations;/**< the maximum iterations */
    int m_Iterations;/**< the current iteration number */
    double m_RAbsTol;/**< the absolution error of the residual */
    double m_RRelTol;/**< the relative error of the residual */
    bool m_IsConverged;/**< converge status */
    double m_Rnorm;/**< the norm of the residual */
    double m_Rnorm0;/**< the initial norm of the residual */
    double m_dUnorm;/**< the norm of the delta u */
    double m_dUnorm0;/**< the initial norm of the delta u */

    Timer m_Timer;/**< timer */

};
