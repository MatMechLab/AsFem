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

    if(!m_inputSystem.isReadOnly()){
        if(m_jobblock.m_jobtype==FEJobType::STATIC){
            runStaticAnalysis();
        }
        else if(m_jobblock.m_jobtype==FEJobType::TRANSIENT){
            runTransientAnalysis();
        }
    }
    else{
        MessagePrinter::printNormalTxt("AsFem has been executed in 'read-only' mode");
        m_timer.endTimer();
        m_timer.printElapseTime("'Simulation' is done");
    }
}