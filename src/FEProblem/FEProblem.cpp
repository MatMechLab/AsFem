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
//+++ Purpose: the FEProblem class of AsFem, the top level of the
//+++          whole program
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

FEProblem::FEProblem(){
    m_timer.resetTimer();
}

void FEProblem::initFEProblem(int args,char *argv[]){
    //***************************************
    // for input file reading
    //***************************************
    m_timer.startTimer();
    m_inputSystem.init(args,argv);

    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Start to read the input file");
    m_inputSystem.readInputFile(m_mesh,m_dofhandler,m_elmtsystem,m_fe,
                                m_bcsystem,m_icsystem,
                                m_projsystem,
                                m_nlsolver,
                                m_timestepping,
                                m_output,
                                m_postprocessor,
                                m_jobblock);
    m_timer.endTimer();
    m_timer.printElapseTime("Input file reading is done",false);
    m_mesh.printBulkMeshInfo();

    //***************************************
    // for dofs init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to create dofs map ...");
    m_dofhandler.createBulkDofsMap(m_mesh,m_elmtsystem);
    m_timer.endTimer();
    m_timer.printElapseTime("Dofs map generation is done",false);
    m_dofhandler.printBulkDofsInfo();

    //***************************************
    // for elmt system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Element system ...");
    m_elmtsystem.init(m_mesh);
    m_timer.endTimer();
    m_timer.printElapseTime("Element system is initialized",false);

    //***************************************
    // for bc system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the BC system ...");
    m_bcsystem.init(m_dofhandler.getMaxDofsPerNode());
    m_timer.endTimer();
    m_timer.printElapseTime("BC system is initialized",false);

    //***************************************
    // for ic system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the IC system ...");
    m_icsystem.init(m_dofhandler.getMaxDofsPerNode());
    m_timer.endTimer();
    m_timer.printElapseTime("IC system is initialized",false);

    //***************************************
    // for FE space init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the FE space ...");
    m_fe.init(m_mesh);
    m_timer.endTimer();
    m_timer.printElapseTime("FE space is initialized",false);

    //***************************************
    // for FE system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the FE system ...");
    m_fesystem.init(m_mesh,m_dofhandler);
    m_timer.endTimer();
    m_timer.printElapseTime("FE system is initialized",false);

    //***************************************
    // for Equation system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Equation system ...");
    m_equationsystem.init(m_dofhandler);
    MessagePrinter::printNormalTxt("  Start to create Sparsity pattern ...");
    m_equationsystem.createSparsityPattern(m_dofhandler);
    MessagePrinter::printNormalTxt("  Sparsity pattern is ready");
    m_timer.endTimer();
    m_timer.printElapseTime("Equation system is initialized",false);

    //***************************************
    // for Solution system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Solution system ...");
    m_solutionsystem.init(m_dofhandler,m_fe);
    m_timer.endTimer();
    m_timer.printElapseTime("Solution system is initialized",false);

    //***************************************
    // for Projection system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Projection system ...");
    m_projsystem.init(m_mesh,m_dofhandler);
    m_timer.endTimer();
    m_timer.printElapseTime("Projection system is initialized",false);

    //***************************************
    // for Nonlinear solver system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the NL solver ...");
    m_nlsolver.init();
    m_timer.endTimer();
    m_timer.printElapseTime("NL solver is initialized",false);

    //***************************************
    // for Postprocess system init
    //***************************************
    m_timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the postprocessor ...");
    m_postprocessor.init();
    m_timer.endTimer();
    m_timer.printElapseTime("Postprocessor is initialized",false);

    //***************************************
    // for fe control init
    //***************************************
    m_fectrlinfo.init();
    m_fectrlinfo.IsDebug=m_jobblock.m_isdebug;
    m_fectrlinfo.IsDepDebug=m_jobblock.m_isdepdebug;


    //***************************************
    // for print out basic info
    //***************************************
    m_elmtsystem.printElmtSystemInfo();
    m_bcsystem.printBCSystemInfo();
    m_icsystem.printICSystemInfo();
    m_projsystem.printProjectionInfo();
    m_fe.printFEInfo();
    m_nlsolver.printSolverInfo();
    m_output.printInfo();
    m_postprocessor.printInfo();
    if(m_jobblock.m_jobtype==FEJobType::TRANSIENT){
        m_timestepping.printInfo();
    }
    m_jobblock.printJobInfo();
}
//*******************************************
void FEProblem::finalize(){
    m_mesh.releaseMemory();
    m_dofhandler.releaseMemory();
    m_fe.releaseMemory();
    m_elmtsystem.releaseMemory();
    m_bcsystem.releaseMemory();
    m_icsystem.releaseMemory();
    m_equationsystem.releaseMemory();
    m_solutionsystem.releaseMemory();
    m_projsystem.releaseMemory();
    m_nlsolver.releaseMemory();
    m_postprocessor.releaseMemory();
    
}