//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FEProblem/FEProblem.h"

void FEProblem::RunTransientAnalysis(){
    _timestepping.SteppingNew(_mesh,_dofHandler,
                        _bcSystem,_icSystem,
                        _elmtSystem,_mateSystem,
                        _equationSystem,_solution,
                        _fe,_feSystem,
                        _outputSystem,_nonlinearsolver,_feCtrlInfo);
                           
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    
}