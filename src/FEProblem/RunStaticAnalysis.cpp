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
//+++ Date   : 2020.12.27
//+++ Purpose: run the static analysis 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::RunStaticAnalysis(){
    MessagePrinter::PrintNormalTxt("Start to do the static FEM analysis ...");
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    _feCtrlInfo._timesteppingtype=TimeSteppingType::STATIC;
    _icSystem.ApplyIC(_mesh,_dofHandler,_solutionSystem._U);
    VecCopy(_solutionSystem._U,_solutionSystem._Unew);//Copy U -> Unew
    VecCopy(_solutionSystem._U,_solutionSystem._Utemp);
    _feSystem.FormBulkFE(FECalcType::InitMaterial,0.0,0.0,_feCtrlInfo.ctan,_mesh,
            _dofHandler,_fe,_elmtSystem,_mateSystem,_solutionSystem,
            _equationSystem._AMATRIX,_equationSystem._RHS);
    _solutionSystem.UpdateMaterials();
    if(_nonlinearSolver.Solve(_mesh,_dofHandler,_elmtSystem,_mateSystem,
        _bcSystem,
        _solutionSystem,_equationSystem,
        _fe,_feSystem,_feCtrlInfo)){
        if(_rank==0){
            _TimerEnd=chrono::high_resolution_clock::now();
            _Duration=Duration(_TimerStart,_TimerEnd);
        }
        char buff[70];string str;
        snprintf(buff,70,"Static analysis finished! [elapse time=%14.6e s]",_Duration);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        
        VecCopy(_solutionSystem._Unew,_solutionSystem._U);//Copy Unew -> U
        VecCopy(_solutionSystem._Unew,_solutionSystem._Utemp);

        if(_feCtrlInfo.IsProjection){
            _feSystem.FormBulkFE(FECalcType::Projection,_feCtrlInfo.dt,_feCtrlInfo.dt,_feCtrlInfo.ctan,
                _mesh,_dofHandler,_fe,_elmtSystem,_mateSystem,
                _solutionSystem,
                _equationSystem._AMATRIX,_equationSystem._RHS);
        }
        _outputSystem.WriteResultToFile(_mesh,_dofHandler,_solutionSystem);
        _postprocessSystem.RunPostprocess(0.0,_mesh,_dofHandler,_fe,_solutionSystem);
        MessagePrinter::PrintStars();
        MessagePrinter::PrintNormalTxt("Write result to "+_outputSystem.GetOutputFileName());
        MessagePrinter::PrintStars();
    }
    else{
        MessagePrinter::PrintNormalTxt("SNES solver failed for your static analysis, please check either your code or your boundary condition");
        MessagePrinter::AsFem_Exit();
    }
}
