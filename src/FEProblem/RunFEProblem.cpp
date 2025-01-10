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
//+++ Date   : 2022.05.11
//+++ Purpose: execute the finite element analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"


void FEProblem::run(){
    int cpus;
    MPI_Comm_size(PETSC_COMM_WORLD,&cpus);
    MessagePrinter::printDashLine(MessageColor::BLUE);
    if(cpus==1){
        MessagePrinter::printNormalTxt(to_string(cpus)+" CPU is used for the simulation",MessageColor::BLUE);
    }
    else{
        MessagePrinter::printNormalTxt(to_string(cpus)+" CPUs are used for the simulation",MessageColor::BLUE);
    }
    MessagePrinter::printDashLine(MessageColor::BLUE);

    if(!m_InputSystem.isReadOnly()){
        if(m_JobBlock.m_JobType==FEJobType::STATIC){
            runStaticAnalysis();
        }
        else if(m_JobBlock.m_JobType==FEJobType::TRANSIENT){
            runTransientAnalysis();
        }
    }
    else{
        MessagePrinter::printNormalTxt("AsFem has been executed in 'read-only' mode");
        m_Timer.endTimer();
        m_Timer.printElapseTime("'Simulation' is done");
    }
}