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
//+++ Purpose: run the static analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::runStaticAnalysis(){
    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Start the static analysis ...");
    MessagePrinter::printStars();
    m_Timer.startTimer();

    // initialize the material
    m_FESystem.formBulkFE(FECalcType::INITMATERIAL,m_FECtrlInfo.T,m_FECtrlInfo.Dt,m_FECtrlInfo.Ctan,
                        m_FECell,m_DofHandler,m_FE,
                        m_ElmtSystem,m_MateSystem,
                        m_SolnSystem,
                        m_EqSystem.m_AMATRIX,m_EqSystem.m_RHS);
    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("Material properties have been initialized");
    MessagePrinter::printDashLine();

    if(m_NLSolver.solve(m_FECell,m_DofHandler,m_FE,
                        m_ElmtSystem,m_MateSystem,
                        m_FESystem,
                        m_BCSystem,
                        m_SolnSystem,m_EqSystem,
                        m_FECtrlInfo)){
        m_Timer.endTimer();
        m_Timer.printElapseTime("Static analysis is done");
        m_ProjSystem.executeProjection(m_FECell,m_DofHandler,m_ElmtSystem,m_MateSystem,m_FE,m_SolnSystem,m_FECtrlInfo);
        m_Output.saveResults2File(-1,m_FECell,m_DofHandler,m_SolnSystem,m_ProjSystem);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printNormalTxt("Save result to "+m_Output.getOutputFileName(),MessageColor::BLUE);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        if(m_PostProcessor.hasPostprocess()){
            m_PostProcessor.prepareCSVFileHeader();
            m_PostProcessor.executePostprocess(m_FECell,m_DofHandler,m_FE,m_MateSystem,m_ProjSystem,m_SolnSystem);
            m_PostProcessor.savePPSResults2CSVFile(0.0);
            MessagePrinter::printDashLine(MessageColor::BLUE);
            MessagePrinter::printNormalTxt("Save postprocess result to "+m_PostProcessor.getCSVFileName(),MessageColor::BLUE);
            MessagePrinter::printDashLine(MessageColor::BLUE);
        }
        MessagePrinter::printStars();
    }
    else{
        MessagePrinter::printErrorTxt("Static analysis fails, please check either your code or your input file and boundary conditions");
        MessagePrinter::exitAsFem();
    }
}