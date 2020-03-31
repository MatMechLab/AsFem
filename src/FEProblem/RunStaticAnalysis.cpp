//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FEProblem/FEProblem.h"

void FEProblem::RunStaticAnalysis(){
    // _nonlinearsolver.Solve(_mesh,_dofHandler,
    //                        _elmtSystem,_mateSystem,
    //                        _bcSystem,_icSystem,_solution,_equationSystem,
    //                        _fe,
    //                        _feSystem);

    // chrono::high_resolution_clock::time_point mystart,myend;
    // if(_rank==0){
    //     mystart=chrono::high_resolution_clock::now();
    // }

    _feCtrlInfo.timesteppingtype=TimeSteppingType::BackWardEuler;
    _feCtrlInfo.ctan[0]=1.0;
    _feCtrlInfo.ctan[1]=1.0;
    _nonlinearsolver.SSolve(_mesh,_dofHandler,
                           _elmtSystem,_mateSystem,
                           _bcSystem,_icSystem,_solution,_equationSystem,
                           _fe,
                           _feSystem,
                           _feCtrlInfo);
    
    if(_feCtrlInfo.IsProjection){
        _feSystem.FormFE(9,_feCtrlInfo.t,_feCtrlInfo.dt,_feCtrlInfo.ctan,
                        _mesh,_dofHandler,_fe,
                        _elmtSystem,_mateSystem,
                        _solution._Unew,_solution._V,
                        _solution._Hist,_solution._HistOld,
                        _solution._Proj,
                        _equationSystem._AMATRIX,_equationSystem._RHS);
        _outputSystem.WriteResultToVTU(_mesh,_dofHandler,_solution._Unew,_solution.GetProjNumPerNode(),_solution.GetProjNameVec(),_solution._Proj);
    }
    else{
        _outputSystem.WriteResultToVTU(_mesh,_dofHandler,_solution._Unew);
    }
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Write result to [%41s] !!!   ***\n",_outputSystem.GetVTUFileName().c_str());
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");

    // if(_rank==0){
    //     myend=chrono::high_resolution_clock::now();
    // }

    // PetscPrintf(PETSC_COMM_WORLD,"*** Nonlinear solver using time=%14.6e [s] \n",
    // chrono::duration_cast<std::chrono::microseconds>(myend-mystart).count()/1.0e6);

    //************************************************************************************
    //*** this is for performance test(just assemble, dont solve the equations!!!)
    //************************************************************************************
    
    // if(_rank==0){
    //     mystart=chrono::high_resolution_clock::now();
    // }

    // for(int i=1;i<=10;i++){
    //     _feSystem.FormFE(6, _feCtrlInfo.t, _feCtrlInfo.dt, _feCtrlInfo.ctan,
    //                      _mesh, _dofHandler, _fe,
    //                      _elmtSystem, _mateSystem,
    //                      _solution._Unew, _solution._V,
    //                      _solution._Hist, _solution._HistOld,
    //                      _solution._Proj,
    //                      _equationSystem._AMATRIX, _equationSystem._RHS);
    // }
    
    
    // if(_rank==0){
    //     myend=chrono::high_resolution_clock::now();
    // }

    // PetscPrintf(PETSC_COMM_WORLD,"*** System assemble using time=%14.6e [s] \n",
    // chrono::duration_cast<std::chrono::microseconds>(myend-mystart).count()/1.0e6);



}