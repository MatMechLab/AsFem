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
//+++ Date   : 2020.12.30
//+++ Purpose: Implement the transient analysis in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::RunTransientAnalysis(){
    MessagePrinter::PrintNormalTxt("Start to do the transient FEM analysis ...");
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    if(_timestepping.Solve(_mesh,_dofHandler,_elmtSystem,_mateSystem,
        _bcSystem,_icSystem,
        _solutionSystem,_equationSystem,
        _fe,_feSystem,
        _outputSystem,
        _postprocessSystem,
        _feCtrlInfo,
        _nonlinearSolver)){
        if(_rank==0){
            _TimerEnd=chrono::high_resolution_clock::now();
            _Duration=Duration(_TimerStart,_TimerEnd);
        }
        MessagePrinter::PrintDashLine(MessageColor::CYAN);
        char buff[70];string str;
        snprintf(buff,70,"Transient analysis finished! [elapse time=%14.6e s]",_Duration);
        str=string(buff);
        MessagePrinter::PrintNormalTxt(str,MessageColor::CYAN);
        MessagePrinter::PrintDashLine(MessageColor::CYAN);
    }
    else{
        MessagePrinter::PrintNormalTxt("Transient analysis failed, please check either your code or your input file");
        MessagePrinter::AsFem_Exit();
    }
}
