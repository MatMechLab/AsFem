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
//+++ Date   : 2022.05.09
//+++ Purpose: run the transient analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::runTransientAnalysis(){

    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Start the transient analysis ...");
    MessagePrinter::printStars();
    m_timer.startTimer();
    if(m_timestepping.solve(m_mesh,m_dofhandler,m_fe,
                            m_elmtsystem,m_matesystem,m_fesystem,
                            m_bcsystem,m_icsystem,
                            m_solutionsystem,m_equationsystem,m_projsystem,
                            m_fectrlinfo,m_nlsolver,m_output,m_postprocessor)){
        m_timer.endTimer();
        MessagePrinter::printStars();
        m_timer.printElapseTime("Transient analysis is done");
        MessagePrinter::printStars();
    }
    else{
        MessagePrinter::printErrorTxt("Transient analysis fails, please check either your code or your input file and boundary conditions");
        MessagePrinter::exitAsFem();
    }
}