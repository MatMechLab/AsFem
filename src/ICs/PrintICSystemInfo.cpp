//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ICs/ICSystem.h"

void ICSystem::PrintICSystemInfo()const{
    if(GetICBlocksNum()>0){
        PetscPrintf(PETSC_COMM_WORLD,"*** Initial condition system information summary:                     ***\n");
        for(auto it:_ICBlockList){
            it.PrintICBlock();
        }
    }
}