//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "DofHandler/DofHandler.h"

void DofHandler::PrintDofInfo() const{
    PetscPrintf(PETSC_COMM_WORLD,"*** DofHandler information summary:                                   ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"***   Dofs=%8d,ActiveDofs=%8d,DofsPerNode=%2d,DofsPerElmt=%3d***\n",GetDofsNum(),GetActiveDofsNum(),GetDofsNumPerNode(),GetDofsNumPerElmt());
    PetscPrintf(PETSC_COMM_WORLD,"***   Dofs name= ");
    for(auto it:_DofNameList){
        PetscPrintf(PETSC_COMM_WORLD,"%-10s ",it.c_str());
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");

}