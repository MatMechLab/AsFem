//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ICs/ICSystem.h"


void ICSystem::AddICBlock(ICBlock &icblock){
    if(_ICBlockList.size()<1){
        _ICBlockList.push_back(icblock);
        _nICBlocks=int(_ICBlockList.size());
    }
    else{
        bool NotInList=true;
        for(unsigned int i=0;i<_ICBlockList.size();i++){
            if(_ICBlockList[i]._ICBlockName==icblock._ICBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            ICBlock temp;

            temp._ICBlockName=icblock._ICBlockName;

            temp._ICTypeName=icblock._ICTypeName;
            temp._ICType=icblock._ICType;

            temp._DofName=icblock._DofName;

            temp._Params=icblock._Params;
            temp._DomainName=icblock._DomainName;
            
            _ICBlockList.push_back(temp);
            _nICBlocks=int(_ICBlockList.size());
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: duplicated [%35s] in [ics] sub block !!!   ***\n",icblock._ICBlockName.c_str());
            Msg_AsFem_Exit();
        }
    }
}