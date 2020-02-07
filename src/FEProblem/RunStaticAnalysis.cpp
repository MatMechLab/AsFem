//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/walkandthinker/AsFem
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
}