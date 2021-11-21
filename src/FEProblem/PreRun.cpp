//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
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
    MessagePrinter::PrintStars();
    MessagePrinter::PrintNormalTxt("Start to read the input file ...");
    _inputSystem.ReadInputFile(_mesh,_dofHandler,
    _elmtSystem,_mateSystem,_bcSystem,_icSystem,_fe,
    _solutionSystem,_outputSystem,
    _postprocessSystem,
    _nonlinearSolver,_timestepping,
    _feJobBlock);
    MessagePrinter::PrintNormalTxt("Input file reading is done !");
    MessagePrinter::PrintStars();
}
//******************************************
void FEProblem::InitAllComponents(){

    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    char buff[70];
    string str;


    //***************************************************************
    //*** for boundary condition system initializing
    //***************************************************************
    _elmtSystem.InitBulkElmtMateInfo(_mateSystem);// set mate index for each elmt block
    _mateSystem.InitBulkMateSystem();// clean all the materials variables
    _bcSystem.InitBCSystem(_mesh);
    if(!_bcSystem.CheckAppliedBCNameIsValid(_mesh)){
        MessagePrinter::PrintErrorTxt("your boundary condition is not correct, please check the boundary name or your mesh file or your input file");
        MessagePrinter::AsFem_Exit();
    }


    snprintf(buff,70,"Start to creat dof map ...");
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    _dofHandler.CreateBulkMeshDofsMap(_mesh,_bcSystem,_elmtSystem);
    if(_rank==0){
        _TimerEnd=chrono::high_resolution_clock::now();
        _Duration=Duration(_TimerStart,_TimerEnd);
    }
    snprintf(buff,70,"  dof map generated ! [elapsed time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);

    
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
    snprintf(buff,70,"  fe space is initialized ! [elapsed time=%14.6e]",_Duration);
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
                            _mesh.GetBulkMeshBulkElmtsNum(),_mesh.GetBulkMeshNodesNum(),
                            _fe._BulkQPoint.GetQpPointsNum());
    
    if(_rank==0){
        _TimerEnd=chrono::high_resolution_clock::now();
        _Duration=Duration(_TimerStart,_TimerEnd);
    }
    snprintf(buff,70,"  solution is initialized ! [elapsed time=%14.6e]",_Duration);
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
    snprintf(buff,70,"  sparse matrix is initialized ! [elapsed time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);

    //***************************************************************
    //*** for nonlinear solver SNES
    //***************************************************************
    snprintf(buff,70,"Start to initialize nonlinear solver system ...");
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(_rank==0){
        _TimerStart=chrono::high_resolution_clock::now();
    }
    _nonlinearSolver.Init();
    if(_rank==0){
        _TimerEnd=chrono::high_resolution_clock::now();
        _Duration=Duration(_TimerStart,_TimerEnd);
    }
    snprintf(buff,70,"  nonlinear solver is initialized ! [elapsed time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);


    //***************************************************************
    //*** for timestepping solver SNES
    //***************************************************************
    if(_feJobType==FEJobType::TRANSIENT){
        snprintf(buff,70,"Start to initialize time stepping system ...");
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        if(_rank==0){
            _TimerStart=chrono::high_resolution_clock::now();
        }
    }
    
    if(_feJobType==FEJobType::TRANSIENT){
        if(_rank==0){
            _TimerEnd=chrono::high_resolution_clock::now();
            _Duration=Duration(_TimerStart,_TimerEnd);
        }
        snprintf(buff,70,"  time stepping is initialized ! [elapsed time=%14.6e]",_Duration);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
    }


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
    snprintf(buff,70,"  fe system is initialized ! [elapsed time=%14.6e]",_Duration);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);

    _postprocessSystem.InitPPSOutput();
    _postprocessSystem.CheckWhetherPPSIsValid(_mesh);

    MessagePrinter::PrintStars(MessageColor::BLUE);
    MessagePrinter::PrintDashLine(MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("Now all the components are ready, AsFem will start the simulation!",MessageColor::BLUE);
    MessagePrinter::PrintDashLine(MessageColor::BLUE);
    MessagePrinter::PrintStars(MessageColor::BLUE);

    _mesh.PrintMeshInfo();
    _dofHandler.PrintAllDofInfo();

    _elmtSystem.PrintElmtSystemInfo();
    _mateSystem.PrintMateSystemInfo();

    _fe.PrintFEInfo();

    _bcSystem.PrintBCSystemInfo();
    _icSystem.PrintICSystemInfo();

    if(_solutionSystem.IsProjection()){
        _solutionSystem.PrintProjectionInfo();
    }

    _outputSystem.PrintInfo();

    _postprocessSystem.PrintPostprocessInfo();

    _nonlinearSolver.PrintInfo();


    _feJobType=_feJobBlock._jobType;
    _feCtrlInfo.IsDebug=_feJobBlock._IsDebug;
    _feCtrlInfo.IsDepDebug=_feJobBlock._IsDepDebug;
    _feCtrlInfo.IsProjection=_solutionSystem.IsProjection();

    if(_feJobType==FEJobType::TRANSIENT){
        _timestepping.PrintTimeSteppingInfo();
    }

    _feJobBlock.PrintJobInfo();

    MessagePrinter::PrintStars(MessageColor::BLUE);
    if(_size==1){
        snprintf(buff,70,"++++++ %8d CPU will be used for the simulation      ++++++!",_size);
    }
    else{
        snprintf(buff,70,"++++++ %8d CPUs will be used for the simulation      ++++++!",_size);
    }
    str=buff;
    MessagePrinter::PrintNormalTxt(str,MessageColor::BLUE);
    MessagePrinter::PrintStars(MessageColor::BLUE);

}
