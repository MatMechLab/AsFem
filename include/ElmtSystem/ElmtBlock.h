//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_ELMTBLOCK_H
#define ASFEM_ELMTBLOCK_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "petsc.h"

#include "ElmtType.h"
#include "MateSystem/MateType.h"

using namespace std;

//******************************************************************
//*** you element block in your input file should look like:
//*** [blockname]
//***   type=element_name(or uel name,i.e. user1, user2)
//***   dofs=dof1 dof2
//***   mate=mateblockname(optional,not the material type name, instead, the material block's name !!!)
//***   domain=alldomain(optional,default one is 'alldomain',thus all the domain will be applied !!!)
//*** [end]

class ElmtBlock{
public:
    string _ElmtBlockName;
    string _ElmtTypeName;
    ElmtType _ElmtType;
    MateType _MateType;
    int      _MateBlockIndex;
    vector<string> _DofNameList;
    vector<PetscInt> _DofIndexList;
    string _MateBlockName;
    string _DomainName;
    bool _IsPrint=true;

    void Reset(){
        _ElmtBlockName.clear();
        _ElmtTypeName.clear();
        _DofNameList.clear();
        _DofIndexList.clear();
        _ElmtType=ElmtType::NullElmt;
        _MateType=MateType::NullMate;
        _MateBlockIndex=1;
        _MateBlockName.clear();
        _DomainName.clear();
        _IsPrint=true;
    }

    void PrintElmtBlock()const{
        PetscPrintf(PETSC_COMM_WORLD,"*** +Element block information:                                       ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***   element block name = [%40s] ***\n",_ElmtBlockName.c_str());
        if(_DofNameList.size()==1){
            PetscPrintf(PETSC_COMM_WORLD,"***   dof name           = %35s        ***\n",_DofNameList[0].c_str());
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"***   dofs name          =");
            for(auto it:_DofNameList){
                PetscPrintf(PETSC_COMM_WORLD,"%10s ",it.c_str());
            }
            PetscPrintf(PETSC_COMM_WORLD,"\n");
        }
        PetscPrintf(PETSC_COMM_WORLD,"***   mate name          = %35s        ***\n",_MateBlockName.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"***   domain name        = %35s        ***\n",_DomainName.c_str());
    }

};


#endif //ASFEM_ELMTBLOCK_H