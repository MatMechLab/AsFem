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
//+++ Purpose: the FEProblem class of AsFem, the top level of the
//+++          whole program
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"
#include "MPIUtils/MPIDataBus.h"

FEProblem::FEProblem(){
    m_Timer.resetTimer();
}

void FEProblem::initFEProblem(int args,char *argv[]){
    //***************************************
    // for input file reading
    //***************************************
    m_Timer.startTimer();
    m_InputSystem.init(args,argv);

    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Start to read the input file");
    m_InputSystem.readInputFile(m_FECell,m_DofHandler,m_ElmtSystem,m_FE,
                                m_BCSystem,m_ICSystem,
                                m_ProjSystem,
                                m_LinearSolver,
                                m_NLSolver,
                                m_TimeStepping,
                                m_Output,
                                m_PostProcessor,
                                m_JobBlock);
    m_Timer.endTimer();
    MessagePrinter::printStars();
    m_Timer.printElapseTime("Input file reading is done",false);

    //***************************************
    // for mesh init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to distribute fecell ...");
    m_FECell.distributeMesh();
    m_Timer.endTimer();
    m_Timer.printElapseTime("FEcell distribution is done",false);

    /**
     * Save the partition info to vtu file, which is named by the input file name
     */
    string FileName;
    FileName = m_InputSystem.getInputFileName();
    m_FECell.saveFECellPartionInfo2VTUFile(FileName.substr(0,FileName.size()-5)+"-partition.vtu");
    MessagePrinter::printDashLine(MessageColor::BLUE);
    MessagePrinter::printNormalTxt("Save partition info to:"+FileName.substr(0,FileName.size()-5)+"-partition.vtu",MessageColor::BLUE);
    MessagePrinter::printDashLine(MessageColor::BLUE);

    if(m_InputSystem.isReadOnly()) return;

    //***************************************
    // for dofs init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to create dofs map ...");
    m_DofHandler.createBulkDofsMap(m_FECell,m_ElmtSystem);
    m_Timer.endTimer();
    m_Timer.printElapseTime("Dofs map generation is done",false);

    //***************************************
    // for elmt system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Element system ...");
    m_ElmtSystem.init(m_FECell);
    m_Timer.endTimer();
    m_Timer.printElapseTime("Element system is initialized",false);


    //***************************************
    // for bc system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the BC system ...");
    m_BCSystem.init(m_DofHandler.getMaxDofsPerNode());
    m_Timer.endTimer();
    m_Timer.printElapseTime("BC system is initialized",false);

    //***************************************
    // for ic system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the IC system ...");
    m_ICSystem.init(m_DofHandler.getMaxDofsPerNode());
    m_Timer.endTimer();
    m_Timer.printElapseTime("IC system is initialized",false);


    //***************************************
    // for FE space init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the FE space ...");
    m_FE.init(m_FECell);
    m_Timer.endTimer();
    m_Timer.printElapseTime("FE space is initialized",false);

    //***************************************
    // for FE system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the FE system ...");
    m_FESystem.init(m_FECell,m_DofHandler);
    m_Timer.endTimer();
    m_Timer.printElapseTime("FE system is initialized",false);

    //***************************************
    // for Equation system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Equation system ...");
    m_EqSystem.init(m_DofHandler);
    MessagePrinter::printNormalTxt("  Start to create Sparsity pattern ...");
    m_EqSystem.createSparsityPattern(m_DofHandler);
    MessagePrinter::printNormalTxt("  Sparsity pattern is ready");
    m_Timer.endTimer();
    m_Timer.printElapseTime("Equation system is initialized",false);

    //***************************************
    // for Solution system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Solution system ...");
    m_SolnSystem.init(m_DofHandler,m_FE);
    m_Timer.endTimer();
    m_Timer.printElapseTime("Solution system is initialized",false);

    //***************************************
    // for Projection system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the Projection system ...");
    m_ProjSystem.init(m_FECell,m_DofHandler);
    m_Timer.endTimer();
    m_Timer.printElapseTime("Projection system is initialized",false);

    //***************************************
    // for Linear solver system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the linear solver ...");
    m_LinearSolver.init();
    m_Timer.endTimer();
    m_Timer.printElapseTime("Linear solver is initialized",false);

    //***************************************
    // for Nonlinear solver system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the NL solver ...");
    m_NLSolver.init(m_LinearSolver);
    m_Timer.endTimer();
    m_Timer.printElapseTime("NL solver is initialized",false);

    //***************************************
    // for Postprocess system init
    //***************************************
    m_Timer.startTimer();
    MessagePrinter::printNormalTxt("Start to initialize the postprocessor ...");
    m_PostProcessor.init();
    m_Timer.endTimer();
    m_Timer.printElapseTime("Postprocessor is initialized",false);

    //***************************************
    // for fe control init
    //***************************************
    m_FECtrlInfo.init();
    m_FECtrlInfo.IsDebug=m_JobBlock.m_IsDebug;
    m_FECtrlInfo.IsDepDebug=m_JobBlock.m_IsDepDebug;


    //***************************************
    // for print out basic info
    //***************************************
    m_FECell.printSummaryInfo();
    m_DofHandler.printBulkDofsInfo();
    m_ElmtSystem.printElmtSystemInfo();
    m_BCSystem.printBCSystemInfo();
    m_ICSystem.printICSystemInfo();
    m_ProjSystem.printProjectionInfo();
    m_FE.printFEInfo();
    m_LinearSolver.printKSPSolverInfo();
    m_NLSolver.printSolverInfo();
    m_Output.printInfo();
    m_PostProcessor.printInfo();
    if(m_JobBlock.m_JobType==FEJobType::TRANSIENT){
        m_TimeStepping.printInfo();
    }
    m_JobBlock.printJobInfo();
}
//*******************************************
void FEProblem::finalize(){
    m_FECell.releaseMemory();
    m_DofHandler.releaseMemory();
    m_FE.releaseMemory();
    m_ElmtSystem.releaseMemory();
    m_BCSystem.releaseMemory();
    m_ICSystem.releaseMemory();
    m_EqSystem.releaseMemory();
    m_SolnSystem.releaseMemory();
    m_ProjSystem.releaseMemory();
    m_LinearSolver.releaseMemory();
    m_NLSolver.releaseMemory();
    m_PostProcessor.releaseMemory();
    
}