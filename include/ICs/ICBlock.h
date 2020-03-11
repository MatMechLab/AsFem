//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_ICBLOCK_H
#define ASFEM_ICBLOCK_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"


#include "ICType.h"

using namespace std;


//**************************************************
//*** format of your bcblock should looks like:
//*** [blockname]
//***   type=ictpye[const,or random, or circle, or rectangle]
//***   dof=dofname
//***   params=
//***   domain=domain_name
//*** [end]
//***************************************************

class ICBlock{
public:
    string _ICBlockName;
    string _ICTypeName;
    ICType _ICType;
    string _DofName;
    PetscInt _DofIndex;
    vector<PetscReal> _Params;
    string _DomainName;// support multiple boundary name list
    void Reset(){
        _ICBlockName.clear();

        _ICType=ICType::NullIC;
        _ICTypeName.clear();

        _DofName.clear();
        _DofIndex=1;

        _Params.clear();

        _DomainName.clear();
    }

    void PrintICBlock()const{
        PetscPrintf(PETSC_COMM_WORLD,"*** +InitialCondition block information:                              ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***   ic block name = [%45s] ***\n",_ICBlockName.c_str());
        
        PetscPrintf(PETSC_COMM_WORLD,"***   ic type name  = %21s                           ***\n",_ICTypeName.c_str());
        
        PetscPrintf(PETSC_COMM_WORLD,"***   dof name      = %-15s, dof index=%2d                   ***\n",_DofName.c_str(),_DofIndex);
        PetscPrintf(PETSC_COMM_WORLD,"***   ic parameters = ");
        for(auto it:_Params){
            PetscPrintf(PETSC_COMM_WORLD,"%-13.5e ",it);
        }
        PetscPrintf(PETSC_COMM_WORLD,"\n");
        PetscPrintf(PETSC_COMM_WORLD,"***   domain name   = %45s   ***\n",_DomainName.c_str());
    }
};


#endif // ASFEM_ICBLOCK_H