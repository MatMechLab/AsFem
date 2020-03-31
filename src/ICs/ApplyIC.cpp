//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ICs/ICSystem.h"

void ICSystem::ApplyIC(Mesh &mesh,DofHandler &dofHandler,Vec &U){
    int iblock,DofIndex;
    ICBlock icblock;
    if(GetICBlocksNum()<1) return;
    for(iblock=1;iblock<=GetICBlocksNum();++iblock){
        DofIndex=dofHandler.GetDofIndexViaName(_ICBlockList[iblock-1]._DofName);
        switch(_ICBlockList[iblock-1]._ICType){
            case ICType::ConstIC:
                ApplyConstIC(_ICBlockList[iblock-1]._Params,DofIndex,mesh,dofHandler,U);
                break;
            case ICType::RandomIC:
                ApplyRandomIC(_ICBlockList[iblock-1]._Params,DofIndex,mesh,dofHandler,U);
                break;
            case ICType::CircleIC:
                ApplyCircleIC(_ICBlockList[iblock-1]._Params,DofIndex,mesh,dofHandler,U);
                break;
            default:
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported initial condition type                   !!!   ***\n");
                Msg_AsFem_Exit();
                break;
        }
    }
    VecAssemblyBegin(U);
    VecAssemblyEnd(U);
}