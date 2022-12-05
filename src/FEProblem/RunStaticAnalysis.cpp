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
//+++ Purpose: run the static analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::runStaticAnalysis(){
    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Start the static analysis ...");
    MessagePrinter::printStars();
    m_timer.startTimer();

    // initialize the material
    m_fesystem.formBulkFE(FECalcType::INITMATERIAL,m_fectrlinfo.t,m_fectrlinfo.dt,m_fectrlinfo.ctan,
                        m_mesh,m_dofhandler,m_fe,
                        m_elmtsystem,m_matesystem,
                        m_solutionsystem,
                        m_equationsystem.m_amatrix,m_equationsystem.m_rhs);
    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("Material properties have been initialized");
    MessagePrinter::printDashLine();

    if(m_nlsolver.solve(m_mesh,m_dofhandler,m_fe,
                        m_elmtsystem,m_matesystem,
                        m_fesystem,
                        m_bcsystem,
                        m_solutionsystem,m_equationsystem,
                        m_fectrlinfo)){

        m_timer.endTimer();
        m_timer.printElapseTime("Static analysis is done");
        m_projsystem.executeProjection(m_mesh,m_dofhandler,m_elmtsystem,m_matesystem,m_fe,m_solutionsystem,m_fectrlinfo);
        m_output.saveResults2File(-1,m_mesh,m_dofhandler,m_solutionsystem,m_projsystem);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printNormalTxt("Save result to "+m_output.getOutputFileName(),MessageColor::BLUE);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        if(m_postprocessor.hasPostprocess()){
            m_postprocessor.prepareCSVFileHeader();
            m_postprocessor.executePostprocess(m_mesh,m_dofhandler,m_fe,m_matesystem,m_projsystem,m_solutionsystem);
            m_postprocessor.savePPSResults2CSVFile(0.0);
            MessagePrinter::printDashLine(MessageColor::BLUE);
            MessagePrinter::printNormalTxt("Save postprocess result to "+m_postprocessor.getCSVFileName(),MessageColor::BLUE);
            MessagePrinter::printDashLine(MessageColor::BLUE);
        }
        MessagePrinter::printStars();
    }
    else{
        MessagePrinter::printErrorTxt("Static analysis fails, please check either your code or your input file and boundary conditions");
        MessagePrinter::exitAsFem();
    }
}