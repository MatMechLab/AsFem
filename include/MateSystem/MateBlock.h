//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_MATEBLOCK_H
#define ASFEM_MATEBLOCK_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"

#include "MateType.h"

using namespace std;

//******************************************************************
//*** you material block in your input file should look like:
//*** [blockname]
//***   type=material_name(or umate name,i.e. user1, user2)
//***   params=
//*** [end]

class MateBlock{
public:
    string _MateBlockName;

    string _MateTypeName;
    MateType _MateType;

    vector<PetscReal> _Params;

    void Reset(){
        _MateBlockName.clear();

        _MateTypeName.clear();
        _MateType=MateType::NullMate;
        _MateBlockName.clear();
        _Params.clear();
    }

    void PrintElmtBlock()const{
        PetscPrintf(PETSC_COMM_WORLD,"*** +Material block information:                                      ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***   material block name = [%40s]***\n",_MateBlockName.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"***   material type       = %35s       ***\n",_MateTypeName.c_str());
        if(_Params.size()==1){
            PetscPrintf(PETSC_COMM_WORLD,"***   parameters          = %13.5e                        ***\n",_Params[0]);
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"***   parameters          =");
            int count=0,j=0;
            for(auto it:_Params){
                PetscPrintf(PETSC_COMM_WORLD,"%13.5e ",it);
                count+=1;
                j+=1;
                if(count%3==0){
                    if(j<int(_Params.size())){
                        PetscPrintf(PETSC_COMM_WORLD," ***\n***                        ");
                    }
                    else{
                        PetscPrintf(PETSC_COMM_WORLD," ***");
                    }
                    count=0;
                }
            }
            if(_Params.size()%3==1){
                PetscPrintf(PETSC_COMM_WORLD,"                             ***");
            }
            else if(_Params.size()%3==2){
                PetscPrintf(PETSC_COMM_WORLD,"               ***");
            }
            PetscPrintf(PETSC_COMM_WORLD,"\n");
        }
    }


};


#endif // ASFEM_MATEBLOCK_H