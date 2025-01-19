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

#include "NonlinearSolver/NewtonRaphsonSolver.h"
#include "TimeStepping/TimeSteppingTool.h"

NewtonRaphsonSolver::NewtonRaphsonSolver() {
    m_MaxIterations=50;
    m_Iterations=0;
    m_RAbsTol=5.0e-7;
    m_RRelTol=1.0e-12;
    m_IsConverged=false;
    m_Rnorm=1.0;
    m_Rnorm0=1.0;
    m_dUnorm=1.0;
    m_dUnorm0=1.0;
}

bool NewtonRaphsonSolver::solve(FECell &t_FECell,
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
    char buff[68];
    m_Timer.startTimer();
    t_SolnSystem.m_Utemp.copyFrom(t_SolnSystem.m_Ucurrent);
    t_SolnSystem.m_Ucopy.copyFrom(t_SolnSystem.m_Ucurrent);
    m_IsConverged=false;
    m_Iterations=0;
    while (m_Iterations<m_MaxIterations && !m_IsConverged) {
        t_BCSystem.applyPresetBoundaryConditions(FECalcType::UPDATEU,
                                             t_FECtrlInfo.Dt+t_FECtrlInfo.T,
                                             t_FECell,
                                             t_DofHandler,
                                             t_SolnSystem.m_Utemp,
                                             t_SolnSystem.m_Ucopy,
                                             t_SolnSystem.m_Uold,
                                             t_SolnSystem.m_Uolder,
                                             t_SolnSystem.m_V,
                                             t_EqSystem.m_AMATRIX,
                                             t_EqSystem.m_RHS);

        computeTimeDerivatives(t_FECtrlInfo,t_SolnSystem);// calculate V and A, based on Utemp

        t_FESystem.formBulkFE(FECalcType::COMPUTERESIDUALANDJACOBIAN,
                              t_FECtrlInfo.Dt+t_FECtrlInfo.T,
                              t_FECtrlInfo.Dt,
                              t_FECtrlInfo.Ctan,
                              t_FECell,
                              t_DofHandler,
                              t_FE,
                              t_ElmtSystem,
                              t_MateSystem,
                              t_SolnSystem,
                              t_EqSystem.m_AMATRIX,
                              t_EqSystem.m_RHS);

        // t_BCSystem.setDirichletPenalty(t_FESystem.getMaxCoefOfKMatrix()*1.0e10);
        t_BCSystem.setDirichletPenalty(1.0e23);
        t_BCSystem.applyBoundaryConditions(FECalcType::COMPUTERESIDUALANDJACOBIAN,
                                           t_FECtrlInfo.Dt+t_FECtrlInfo.T,
                                           t_FECtrlInfo.Ctan,
                                           t_FECell,
                                           t_DofHandler,
                                           t_FE,
                                           t_SolnSystem.m_Utemp,
                                           t_SolnSystem.m_Ucopy,
                                           t_SolnSystem.m_Uold,
                                           t_SolnSystem.m_Uolder,
                                           t_SolnSystem.m_V,
                                           t_EqSystem.m_AMATRIX,
                                           t_EqSystem.m_RHS);
        t_EqSystem.m_AMATRIX*=-1.0;// in this NR iteration, the K=-dR/dU, which is different from the one in SNES solver !!!
        t_LinearSolver.solve(t_EqSystem.m_AMATRIX,t_EqSystem.m_RHS,t_SolnSystem.m_dU);
        t_SolnSystem.m_Utemp+=t_SolnSystem.m_dU;
        m_Iterations+=1;
        m_Rnorm=t_EqSystem.m_RHS.getNorm();
        m_dUnorm=t_SolnSystem.m_dU.getNorm();
        if (m_Iterations==1) {
            m_Rnorm0=m_Rnorm;
            m_dUnorm0=m_dUnorm;
        }
        if (t_FECtrlInfo.IsDepDebug) {
            snprintf(buff,68,"  AsFem solver: iters=%4d, |R|=%12.5e, |dU|=%12.5e",m_Iterations,m_Rnorm,m_dUnorm);
            MessagePrinter::printNormalTxt(buff);
        }
        if (m_Rnorm<m_RAbsTol||m_Rnorm<m_RRelTol*m_Rnorm0) {
            m_IsConverged=true;
            break;
        }
    }
    if (!t_FECtrlInfo.IsDepDebug) {
        snprintf(buff,68,"  AsFem solver: iters=%4d, |R0|=%12.5e, |R|=%12.5e",m_Iterations,m_Rnorm0,m_Rnorm);
        MessagePrinter::printNormalTxt(buff);
    }
    t_SolnSystem.m_Ucurrent.copyFrom(t_SolnSystem.m_Utemp);
    m_Timer.endTimer();
    m_Timer.printElapseTime("AsFem solver is done");
    return m_IsConverged;
}

void NewtonRaphsonSolver::printSolverInfo()const {
    MessagePrinter::printNormalTxt("Nonlinear (ASFEM) solver information summary:");
    char buff[70];
    string str;


    str="  Solver type= classical newton-raphson";
    MessagePrinter::printNormalTxt(str);

    snprintf(buff,70,"  Max iterations=%3d,",m_MaxIterations);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    snprintf(buff,70,"  Abusolute |R| tolerance=%14.5e",m_RAbsTol);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    snprintf(buff,70,"  Relative |R| tolerance=%14.5e",m_RRelTol);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    MessagePrinter::printStars();
}

