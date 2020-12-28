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
//+++ Purpose: read the input file before we run FEM analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

void FEProblem::ReadInputFile(){
    _inputSystem.ReadInputFile(_mesh,_dofHandler,
    _elmtSystem,_mateSystem,_bcSystem,_icSystem,_fe,
    _solutionSystem,_outputSystem,_nonlinearSolver);

    _mesh.PrintMeshInfo();
    _dofHandler.PrintAllDofInfo();

    _elmtSystem.InitBulkElmtMateInfo(_mateSystem);
    _elmtSystem.PrintElmtSystemInfo();

    _mateSystem.InitBulkMateSystem();// clean all the materials variables
    _mateSystem.PrintMateSystemInfo();

    _fe.PrintFEInfo();

    _bcSystem.PrintBCSystemInfo();
}
//******************************************
void FEProblem::InitAllComponents(){

    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    char buff[70];
    string str;

    snprintf(buff,70,"Start to creat dof map ...");
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    _dofHandler.CreateBulkDofsMap(_mesh,_bcSystem,_elmtSystem);
    if(_rank==0){
        _TimerEnd=chrono::high_resolution_clock::now();
        _Duration=Duration(_TimerStart,_TimerEnd);
    }
    snprintf(buff,70,"  dof map generated ! [time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);

    //***************************************************************
    //*** for boundary condition system initializing
    //***************************************************************
    _bcSystem.InitBCSystem(_mesh);
    
    //***************************************************************
    //*** for FE space initializing
    //***************************************************************
    snprintf(buff,70,"Start to initialize FE space ...");
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    _fe.InitFE(_mesh);
    if(_rank==0){
        _TimerEnd=chrono::high_resolution_clock::now();
        _Duration=Duration(_TimerStart,_TimerEnd);
    }
    snprintf(buff,70,"  fe space is initialized ! [time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    

    //***************************************************************
    //*** for solution system initializing
    //***************************************************************
    snprintf(buff,70,"Start to initialize solution system ...");
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    _solutionSystem.SetHistNumPerGPoint(10);
    _solutionSystem.InitSolution(_dofHandler.GetActiveDofsNum(),
                            _mesh.GetBulkMeshNodesNum(),_mesh.GetBulkMeshBulkElmtsNum(),
                            _fe._BulkQPoint.GetQpPointsNum());
    
    if(_rank==0){
        _TimerEnd=chrono::high_resolution_clock::now();
        _Duration=Duration(_TimerStart,_TimerEnd);
    }
    snprintf(buff,70,"  solution is initialized ! [time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);

    //***************************************************************
    //*** for equation system initializing
    //***************************************************************
    snprintf(buff,70,"Start to initialize sparse matrix system ...");
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    _equationSystem.InitEquationSystem(_dofHandler.GetActiveDofsNum(),_dofHandler.GetMaxRowNNZ());
    _equationSystem.CreateSparsityPattern(_dofHandler);
    if(_rank==0){
        _TimerEnd=chrono::high_resolution_clock::now();
        _Duration=Duration(_TimerStart,_TimerEnd);
    }
    snprintf(buff,70,"  sparse matrix is initialized ! [time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);

    //***************************************************************
    //*** for nonlinear solver
    //*** it is already initialized in input system !!!
    //***************************************************************


    //***************************************************************
    //*** for FE system initializing
    //***************************************************************
    snprintf(buff,70,"Start to initialize FE system ...");
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    _feSystem.InitBulkFESystem(_mesh,_dofHandler,_fe,_solutionSystem);
    if(_rank==0){
        _TimerEnd=chrono::high_resolution_clock::now();
        _Duration=Duration(_TimerStart,_TimerEnd);
    }
    snprintf(buff,70,"  fe system is initialized ! [time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);

    MessagePrinter::PrintStars();

}