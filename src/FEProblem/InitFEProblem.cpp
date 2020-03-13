//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FEProblem/FEProblem.h"

void FEProblem::InitFEProblem(){

    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    if(_rank==0){
        _TimerStartOfFEProInit=chrono::high_resolution_clock::now();
    }
    
    PetscPrintf(PETSC_COMM_WORLD,"*** Start to initialize the FE problem ...                            ***\n");
    
    //***********************************************************
    //*** init the mesh system, set their related properties
    //***********************************************************
    PetscPrintf(PETSC_COMM_WORLD,"***   start to initialize the mesh properties ...                     ***\n");
    
    if(_rank==0) _TimerStartOfElmtInit=chrono::high_resolution_clock::now();
    
    _mesh.SetElmtInfo(_elmtSystem);
    
    if(_rank==0){
         _TimerEndOfElmtInit=chrono::high_resolution_clock::now();
         _DurationOfElmtInit=Duration(_TimerStartOfElmtInit,_TimerEndOfElmtInit);
    }
    
    PetscPrintf(PETSC_COMM_WORLD,"***   mesh properties initialized             !!!===>[%13.5e s]***\n",_DurationOfElmtInit);
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    
    //*************************************************
    // initialize the dof system
    //*************************************************
    if(_rank==0){
        _TimerStartOfDofInit=chrono::high_resolution_clock::now();
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   start to initialize the dof handler ...                         ***\n");
    _dofHandler.CreateDofMap(_mesh,_bcSystem);

    if(_rank==0){
        _TimerEndOfDofInit=chrono::high_resolution_clock::now();
        _DurationOfDofInit=Duration(_TimerStartOfDofInit,_TimerEndOfDofInit);
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   dof handler initialized                 !!!===>[%13.5e s]***\n",_DurationOfDofInit);
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    
    
    
    //*************************************************
    //*** Init sparse matrix system
    //*************************************************
    if(_rank==0){
        _TimerStartOfSparseInit=chrono::high_resolution_clock::now();
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   start to initialize the sparse equation system ...              ***\n");
    _equationSystem.InitEquationSystem(_dofHandler.GetActiveDofsNum(),_dofHandler.GetMaxRowNNZ());
    _equationSystem.CreateSparsityPattern(_dofHandler);
    
    if(_rank==0){
        _TimerEndOfSparseInit=chrono::high_resolution_clock::now();
        _DurationOfSparseInit=Duration(_TimerStartOfSparseInit,_TimerEndOfSparseInit);
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   sparse equation system initialized      !!!===>[%13.5e s]***\n",_DurationOfSparseInit);
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    
    //*************************************************
    //*** Init FE space
    //*************************************************
    if(_rank==0){
        _TimerStartOfFESysInit=chrono::high_resolution_clock::now();
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   start to initialize the FE space system ...                     ***\n");
    _equationSystem.InitEquationSystem(_dofHandler.GetActiveDofsNum(),_dofHandler.GetMaxRowNNZ());
    _equationSystem.CreateSparsityPattern(_dofHandler);
    _fe.InitFE(_mesh);
    
    if(_rank==0){
        _TimerEndOfFESysInit=chrono::high_resolution_clock::now();
        _DurationOfFESysInit=Duration(_TimerStartOfFESysInit,_TimerEndOfFESysInit);
    }

    PetscPrintf(PETSC_COMM_WORLD,"***   FE space system initialized             !!!===>[%13.5e s]***\n",_DurationOfFESysInit);
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    
    //*************************************************
    //*** Init solution system
    //*************************************************
    if(_rank==0){
        _TimerStartOfSolInit=chrono::high_resolution_clock::now();
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   start to initialize the solution system ...                     ***\n");
    
    _solution.AddSolutionNameFromVec(_dofHandler.GetDofNameList());
    _solution.InitSolution(_dofHandler.GetActiveDofsNum(),
                           _mesh.GetBulkElmtsNum(),
                           _mesh.GetNodesNum(),
                           _fe._qp_bulk.GetQpPointsNum());
    if(_rank==0){
        _TimerEndOfSolInit=chrono::high_resolution_clock::now();
        _DurationOfSolInit=Duration(_TimerStartOfSolInit,_TimerEndOfSolInit);
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   solution system initialized             !!!===>[%13.5e s]***\n",_DurationOfSolInit);
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    
    
    //*************************************************
    //*** Init FE System
    //*************************************************
    if(_rank==0){
        _TimerStartOfFESysInit=chrono::high_resolution_clock::now();
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   start to initialize the FE system ...                           ***\n");
    _feSystem.InitFESystem(_mesh,_dofHandler,_fe,_solution);
    _mateSystem.InitMateSystem();
    _bcSystem.InitBCSystem(_mesh);
    
    if(_rank==0){
        _TimerEndOfFESysInit=chrono::high_resolution_clock::now();
        _DurationOfFESysInit=Duration(_TimerStartOfFESysInit,_TimerEndOfFESysInit);
    }

    PetscPrintf(PETSC_COMM_WORLD,"***   FE system initialized                   !!!===>[%13.5e s]***\n",_DurationOfFESysInit);
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    
    //*************************************************
    //*** Init Nonlinear solver system
    //*************************************************
    if(_rank==0){
        _TimerStartOfSolverInit=chrono::high_resolution_clock::now();
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   start to initialize the Nonlinear solver system ...             ***\n");
    _nonlinearsolver.Init(_nonlinearsolverblock);
    _timesteppingblock._interval=_jobBlock._Interval;
    _timestepping.Init(_timesteppingblock);
    _feCtrlInfo.timesteppingtype=_timesteppingblock._TimeSteppingMethod;
    // _timestepping.InitSolver(_nonlinearsolverblock);
    
    if(_rank==0){
        _TimerEndOfSolverInit=chrono::high_resolution_clock::now();
        _DurationOfSolverInit=Duration(_TimerStartOfFESysInit,_TimerEndOfFESysInit);
    }

    PetscPrintf(PETSC_COMM_WORLD,"***   Nonlinear solver system initialized     !!!===>[%13.5e s]***\n",_DurationOfSolverInit);
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");


    //*************************************************
    //*** Init output system
    //*************************************************
    if(_rank==0){
        _TimerStartOfOutput=chrono::high_resolution_clock::now();
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   start to initialize the output system ...                       ***\n");
    
    // _outputSystem.SetInputFileName(_inputSystem.GetInputFileName());
    _outputSystem.InitOutputStream();
    
    if(_rank==0){
        _TimerEndOfOutput=chrono::high_resolution_clock::now();
        _DurationOfOutput=Duration(_TimerStartOfFESysInit,_TimerEndOfFESysInit);
    }

    PetscPrintf(PETSC_COMM_WORLD,"***   Output system initialized               !!!===>[%13.5e s]***\n",_DurationOfOutput);
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");

    //*************************************************
    //*** Finish all the initialization
    //*************************************************
    if(_rank==0){
        _TimerEndOfFEProInit=chrono::high_resolution_clock::now();
        _DurationOfFEProInit=Duration(_TimerStartOfFEProInit,_TimerEndOfFEProInit);
    }
    PetscPrintf(PETSC_COMM_WORLD,"*** FE problem initialized                    !!!===>[%13.5e s]***\n",_DurationOfFEProInit);


    _mesh.PrintMeshInfo();
    _fe.PrintQPointInfo();
    _dofHandler.PrintDofInfo();
    _bcSystem.PrintBCSystemInfo();
    _icSystem.PrintICSystemInfo();
    _elmtSystem.PrintElmtSystemInfo();
    _mateSystem.PrintMateSystemInfo();

    if(_jobBlock._JobType==JobType::StaticJob){
        _nonlinearsolverblock.PrintNonlinearSolverBlock();
        _feCtrlInfo.ctan[0]=1.0;_feCtrlInfo.ctan[1]=1.0;
        _feCtrlInfo.CurrentStep=0;_feCtrlInfo.dt=1.0;
        _feCtrlInfo.t=1.0;
    }
    else if(_jobBlock._JobType==JobType::TransientJob){
        _timesteppingblock.PrintTimeSteppingBlock();
        _nonlinearsolverblock.PrintNonlinearSolverBlock();
    }

    _outputSystem.PrintOutputSystem();
    
    _jobBlock.PrintJobBlockInfo();
    
    _feCtrlInfo.IsDebug=_jobBlock._IsDebug;
    _feCtrlInfo.IsDepDebug=_jobBlock._IsDepDebug;
    _feCtrlInfo.IsProjection=_jobBlock._IsProjection;
    
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");

}