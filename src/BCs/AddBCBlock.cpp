//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "BCs/BCSystem.h"

void BCSystem::AddBCBlock(BCBlock &bcblock){
    if(_BCBlockList.size()<1){
        _BCBlockList.push_back(bcblock);
        _nBCBlocks=int(_BCBlockList.size());
    }
    else{
        bool NotInList=true;
        for(unsigned int i=0;i<_BCBlockList.size();i++){
            if(_BCBlockList[i]._BCBlockName==bcblock._BCBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            BCBlock temp;
            temp._BCBlockName=bcblock._BCBlockName;

            temp._BCTypeName=bcblock._BCTypeName;
            temp._BCType=bcblock._BCType;

            temp._DofName=bcblock._DofName;
            temp._DofIndex=bcblock._DofIndex;

            temp._BCValue=bcblock._BCValue;

            temp._BoundaryNameList=bcblock._BoundaryNameList;

            temp._IsTimeDependent=bcblock._IsTimeDependent;
            _BCBlockList.push_back(temp);
            _nBCBlocks=int(_BCBlockList.size());
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: duplicate [%22s] in [bcs] subblock !!!   ***\n",bcblock._BCBlockName.c_str());
            Msg_AsFem_Exit();
        }
    }
}