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

#include "NonlinearSolver/NonlinearSolver.h"

NonlinearSolver::NonlinearSolver(){
    m_NlSolverBlock.init();
}

void NonlinearSolver::init(LinearSolver &lsolver){
    setFromNonlinearSolverBlock(m_NlSolverBlock);
    if (m_NlSolverBlock.m_NlSolverType==NonlinearSolverType::ASFEMNR) {
        NewtonRaphsonSolver::setMaxIterationNum(m_NlSolverBlock.m_MaxIters);
        NewtonRaphsonSolver::setRAbsTolerance(m_NlSolverBlock.m_AbsTolR);
        NewtonRaphsonSolver::setRRelTolerance(m_NlSolverBlock.m_RelTolR);
    }
    else {
        initSolver(lsolver);
    }
}

int NonlinearSolver::getIterationNum()const {
    if (m_NlSolverBlock.m_NlSolverType==NonlinearSolverType::ASFEMNR) {
        return NewtonRaphsonSolver::getNRIterationNum();
    }
    else {
        return SNESSolver::getIterationNum();
    }
    return -1;
}

bool NonlinearSolver::solve(FECell &t_FECell,
                            DofHandler &t_DofHandler,
                            FE &t_FE,
                            ElmtSystem &t_ElmtSystem,
                            MateSystem &t_MateSystem,
                            FESystem &t_FESystem,
                            BCSystem &t_BCSystem,
                            SolutionSystem &t_SolnSystem,
                            EquationSystem &t_EqSystem,
                            LinearSolver &t_LinearSolver,
                            FEControlInfo &t_FECtrlInfo) {
    if (m_NlSolverBlock.m_NlSolverType==NonlinearSolverType::ASFEMNR) {
        return NewtonRaphsonSolver::solve(t_FECell,
                                   t_DofHandler,
                                   t_FE,
                                   t_ElmtSystem,
                                   t_MateSystem,
                                   t_FESystem,
                                   t_BCSystem,
                                   t_SolnSystem,
                                   t_EqSystem,
                                   t_LinearSolver,
                                   t_FECtrlInfo);
    }
    else {
        return SNESSolver::solve(t_FECell,
                          t_DofHandler,
                          t_FE,
                          t_ElmtSystem,
                          t_MateSystem,
                          t_FESystem,
                          t_BCSystem,
                          t_SolnSystem,
                          t_EqSystem,
                          t_LinearSolver,
                          t_FECtrlInfo);
    }
    return false;
}

void NonlinearSolver::printSolverInfo()const {
    if (m_NlSolverBlock.m_NlSolverType==NonlinearSolverType::ASFEMNR) {
        NewtonRaphsonSolver::printSolverInfo();
    }
    else {
        SNESSolver::printSolverInfo();
    }
}

void NonlinearSolver::releaseMemory() {
    if (m_NlSolverBlock.m_NlSolverType!=NonlinearSolverType::ASFEMNR) {
        SNESSolver::releaseMemory();
    }
}