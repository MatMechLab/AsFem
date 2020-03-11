//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::PrintMateSystemInfo()const{
    if(_nMateBlocks>0){
        PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
        PetscPrintf(PETSC_COMM_WORLD,"*** Material system information summary:                              ***\n");
        for(auto it:_MateBlockList){
            it.PrintElmtBlock();
        }
    }
}