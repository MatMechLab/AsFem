//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.27
//+++ Purpose: run the related FEM analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::Run(){
    ReadInputFile();
    InitAllComponents();

    if(_feJobType==FEJobType::STATIC){
        RunStaticAnalysis();
    }
    else if(_feJobType==FEJobType::TRANSIENT){
        RunTransientAnalysis();
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported FEM job type, please check your input file");
        MessagePrinter::AsFem_Exit();
    }
}