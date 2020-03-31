//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_BCBLOCK_H
#define ASFEM_BCBLOCK_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"


#include "BCType.h"

using namespace std;


//**************************************************
//*** format of your bcblock should looks like:
//*** [blockname]
//***   type=bctpye[dirichlet,or neumann, or userxx]
//***   dof=dofname
//***   value=values or t or 1.0*t or t*1.0
//***   boundary=boundary_name
//*** [end]
//***************************************************

class BCBlock{
public:
    string _BCBlockName;
    string _BCTypeName;
    BCType _BCType;
    string _DofName;
    PetscInt _DofIndex;
    PetscReal _BCValue;         //if it is t or value*t, then only the value is stored,
    bool _IsTimeDependent=false;// TimeDependent will be true
    vector<string> _BoundaryNameList;// support multiple boundary name list
    void Reset(){
        _BCBlockName.clear();
        _BCType=BCType::NullBC;
        _BCTypeName.clear();
        _DofName.clear();
        _DofIndex=1;
        _BCValue=0.0;
        _IsTimeDependent=false;
        _BoundaryNameList.clear();
    }

    void PrintBCBlock()const{
        PetscPrintf(PETSC_COMM_WORLD,"*** +Boundary block information:                                      ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***   boundary block name = [%40s]***\n",_BCBlockName.c_str());
        if(_IsTimeDependent){
            PetscPrintf(PETSC_COMM_WORLD,"***   boundary type name  = %15s (Time dependent)          ***\n",_BCTypeName.c_str());
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"***   boundary type name  = %15s                           ***\n",_BCTypeName.c_str());
        }
        PetscPrintf(PETSC_COMM_WORLD,"***   dof name            = %15s, dof index=%2d             ***\n",_DofName.c_str(),_DofIndex);
        PetscPrintf(PETSC_COMM_WORLD,"***   boundary value      = %14.6e                            ***\n",_BCValue);
        PetscPrintf(PETSC_COMM_WORLD,"***   boundary name       =");
        int count=0;
        for(auto it:_BoundaryNameList){
            PetscPrintf(PETSC_COMM_WORLD,"%-12s ",it.c_str());
            count+=1;
            if(count%5==0&&_BoundaryNameList.size()%5!=0){
                PetscPrintf(PETSC_COMM_WORLD,"\n***                       =");
            }
        }
        PetscPrintf(PETSC_COMM_WORLD,"\n");
    }
};


#endif // ASFEM_BCBLOCK_H