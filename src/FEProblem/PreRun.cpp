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
}
//******************************************
void FEProblem::InitAllComponents(){
    _dofHandler.CreateBulkDofsMap(_mesh,_bcSystem,_elmtSystem);

    _elmtSystem.InitBulkElmtMateInfo(_mateSystem);

    _bcSystem.InitBCSystem(_mesh);
    
    _fe.InitFE(_mesh);
    
    _solutionSystem.SetHistNumPerGPoint(10);
    _solutionSystem.InitSolution(_dofHandler.GetActiveDofsNum(),
                            _mesh.GetBulkMeshNodesNum(),_mesh.GetBulkMeshBulkElmtsNum(),
                            _fe._BulkQPoint.GetQpPointsNum());

    _equationSystem.InitEquationSystem(_dofHandler.GetActiveDofsNum(),_dofHandler.GetMaxRowNNZ());
    _equationSystem.CreateSparsityPattern(_dofHandler);

}