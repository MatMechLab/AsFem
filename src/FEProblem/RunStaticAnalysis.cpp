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
//+++ Purpose: run the static analysis 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::RunStaticAnalysis(){
    MessagePrinter::PrintNormalTxt("Start to do the static FEM analysis ...");
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }

    if(_nonlinearSolver.Solve(_mesh,_dofHandler,_elmtSystem,_mateSystem,
        _bcSystem,_icSystem,
        _solutionSystem,_equationSystem,
        _fe,_feSystem,_feCtrlInfo)){
        if(_rank==0){
            _TimerEnd=chrono::high_resolution_clock::now();
            _Duration=Duration(_TimerStart,_TimerEnd);
        }
        char buff[70];string str;
        snprintf(buff,70,"Static analysis finished! [elapse time=%14.6e]",_Duration);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        if(_feCtrlInfo.IsProjection){
            _feSystem.FormBulkFE(FECalcType::Projection,_feCtrlInfo.dt,_feCtrlInfo.dt,_feCtrlInfo.ctan,
                _mesh,_dofHandler,_fe,_elmtSystem,_mateSystem,
                _solutionSystem._Unew,_solutionSystem._V,
                _solutionSystem._Hist,_solutionSystem._HistOld,_solutionSystem._Proj,
                _equationSystem._AMATRIX,_equationSystem._RHS);

        
            _outputSystem.WriteResultToFile(_mesh,_dofHandler,
            _solutionSystem._Unew,_solutionSystem.GetProjNumPerNode(),_solutionSystem.GetProjNameVec(),_solutionSystem._Proj);
        }
        else{
             _outputSystem.WriteResultToFile(_mesh,_dofHandler,_solutionSystem._Unew);
        }
        MessagePrinter::PrintStars();
        MessagePrinter::PrintNormalTxt("Write result to "+_outputSystem.GetOutputFileName());
        MessagePrinter::PrintStars();
    }
    else{
        MessagePrinter::PrintNormalTxt("SNES solver failed for your static analysis, please check either your code or your boundary condition");
        MessagePrinter::AsFem_Exit();
    }
}