//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "BCs/BCSystem.h"

void BCSystem::PrintBCSystemInfo()const{
    if(_nBCBlocks>0){
        PetscPrintf(PETSC_COMM_WORLD,"*** Boundary condition system information summary:                    ***\n");
        for(auto it:_BCBlockList){
            it.PrintBCBlock();
        }
    }
}