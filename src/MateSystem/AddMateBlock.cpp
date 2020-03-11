//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::AddMateBlock(MateBlock &mateblock){
    if(_MateBlockList.size()<1){
        _MateBlockList.push_back(mateblock);
        _nMateBlocks=int(_MateBlockList.size());
    }
    else{
        bool NotInList=true;
        for(unsigned int i=0;i<_MateBlockList.size();i++){
            if(_MateBlockList[i]._MateBlockName==mateblock._MateBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            // MateBlock temp;
            // temp._MateBlockName=mateblock._MateBlockName;
            // temp._MateTypeName=mateblock._MateTypeName;
            // temp._Params=mateblock._Params;
            
            _MateBlockList.push_back(mateblock);
            _nMateBlocks=int(_MateBlockList.size());
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: duplicated [%35s] in [mate] sub block !!!   ***\n",mateblock._MateBlockName.c_str());
            Msg_AsFem_Exit();
        }
    }
}