//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::AddElmtBlock(ElmtBlock &elmtBlock){
    if(_ElmtBlockList.size()<1){
        _ElmtBlockList.push_back(elmtBlock);
        _nElmtBlocks=int(_ElmtBlockList.size());
    }
    else{
        bool NotInList=true;
        for(unsigned int i=0;i<_ElmtBlockList.size();i++){
            if(_ElmtBlockList[i]._ElmtBlockName==elmtBlock._ElmtBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            ElmtBlock tempblock;
            
            tempblock._ElmtBlockName=elmtBlock._ElmtBlockName;

            tempblock._ElmtTypeName=elmtBlock._ElmtTypeName;
            tempblock._ElmtType=elmtBlock._ElmtType;

            tempblock._DofNameList=elmtBlock._DofNameList;
            tempblock._DofIndexList=elmtBlock._DofIndexList;

            tempblock._DomainName=elmtBlock._DomainName;
            tempblock._MateBlockName=elmtBlock._MateBlockName;

            tempblock._IsPrint=elmtBlock._IsPrint;
            
            _ElmtBlockList.push_back(tempblock);
            _nElmtBlocks=int(_ElmtBlockList.size());
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: duplicated [%35s] in [elmts] block !!!   ***\n",elmtBlock._ElmtBlockName.c_str());
            Msg_AsFem_Exit();
        }
    }
}