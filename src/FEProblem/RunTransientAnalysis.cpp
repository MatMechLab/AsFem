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
//+++ Date   : 2022.05.09
//+++ Purpose: run the transient analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::runTransientAnalysis(){

    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Start the transient analysis ...");
    MessagePrinter::printStars();
    m_Timer.startTimer();
    if(m_TimeStepping.solve(m_FECell,
                            m_DofHandler,
                            m_FE,
                            m_ElmtSystem,
                            m_MateSystem,
                            m_FESystem,
                            m_BCSystem,
                            m_ICSystem,
                            m_SolnSystem,
                            m_EqSystem,
                            m_ProjSystem,
                            m_FECtrlInfo,
                            m_LinearSolver,
                            m_NLSolver,
                            m_Output,
                            m_PostProcessor)){
        m_Timer.endTimer();
        MessagePrinter::printStars();
        m_Timer.printElapseTime("Transient analysis is done");
        MessagePrinter::printStars();
    }
    else{
        MessagePrinter::printErrorTxt("Transient analysis fails, please check either your code or your input file and boundary conditions");
        MessagePrinter::exitAsFem();
    }
}