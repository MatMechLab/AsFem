//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.29
//+++ Purpose: write results to files in different format
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "OutputSystem/OutputSystem.h"

void OutputSystem::WriteResultToFile(const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem){
    if(_OutputType==OutputType::VTU){
        WriteResult2VTU(mesh,dofHandler,solutionSystem);
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported output file format, we will update this in the future");
        MessagePrinter::AsFem_Exit();
    }
}
void OutputSystem::WriteResultToFile(const int &step,const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem){
    if(_OutputType==OutputType::VTU){
        WriteResult2VTU(step,mesh,dofHandler,solutionSystem);
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported output file format, we will update this in the future");
        MessagePrinter::AsFem_Exit();
    }
}
    