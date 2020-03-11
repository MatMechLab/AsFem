//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Mesh/Mesh.h"

void Mesh::PrintMeshInfo() const{
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Mesh information summary:                                         ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"***  nodes=%9d, elmts=%9d, nodesperbulkelmt=%3d           ***\n",GetNodesNum(),
                                                                                           GetElmtsNum(),
                                                                                           GetNodesNumPerBulkElmt());
    PetscPrintf(PETSC_COMM_WORLD,"***  max dim=%2d, min dim=%2d, phygroup=%5d, meshtype=%6s, order=%1d ***\n",GetDim(),
                                                                                             GetMinDim(),
                                                                                             GetPhysicalGroupNum(),
                                                                                             _BulkMeshTypeName.c_str(),
                                                                                             GetMeshOrder());
    PetscPrintf(PETSC_COMM_WORLD,"***  physical id                    phsical Name             elmts    ***\n");
    for(PetscInt i=0;i<GetPhysicalGroupNum();++i){
        PetscPrintf(PETSC_COMM_WORLD,"***   %6d      %30s          %8d    ***\n",_PhysicIDToNameList[i].first,
                                                                        _PhysicIDToNameList[i].second.c_str(),
                                                                        _PhysicNameToElmtIndexSet[i].second.size());
    }
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
}
//****************************************************
void Mesh::PrintMeshDetailInfo() const{
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
    PetscPrintf(PETSC_COMM_WORLD,"*** Mesh information summary:                                         ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"***  nodes=%9d, elmts=%9d, nodesperbulkelmt=%3d           ***\n",GetNodesNum(),
                                                                                           GetElmtsNum(),
                                                                                           GetNodesNumPerBulkElmt());
    PetscPrintf(PETSC_COMM_WORLD,"***  max dim=%2d, min dim=%2d, phygroup=%5d, meshtype=%6s, order=%1d ***\n",GetDim(),
                                                                                             GetMinDim(),
                                                                                             GetPhysicalGroupNum(),
                                                                                             _BulkMeshTypeName.c_str(),
                                                                                             GetMeshOrder());
    PetscPrintf(PETSC_COMM_WORLD,"***  physical id                    phsical name             elmts    ***\n");
    for(PetscInt i=0;i<GetPhysicalGroupNum();++i){
        PetscPrintf(PETSC_COMM_WORLD,"***   %6d      %30s          %8d    ***\n",_PhysicIDToNameList[i].first,
                                                                        _PhysicIDToNameList[i].second.c_str(),
                                                                        _PhysicNameToElmtIndexSet[i].second.size());
    }

    PetscPrintf(PETSC_COMM_WORLD,"***  Physical group information (ID and element ID)                   ***\n");
    for(auto it:_PhysicNameToElmtIndexSet){
        for(PetscInt e=0;e<(PetscInt)it.second.size();++e){
            PetscPrintf(PETSC_COMM_WORLD,"***  phyname=%25s, element id=%9d          ***\n",it.first.c_str(),it.second[e]);
        }
    }  

    PetscPrintf(PETSC_COMM_WORLD,"***  element connectivity information(element id: node index):        ***\n");
    for(PetscInt e=1;e<=GetElmtsNum();++e){
        PetscPrintf(PETSC_COMM_WORLD,"***  elmt id=%9d:",e);
        for(PetscInt i=1;i<=GetIthElmtNodesNum(e);++i){
            if(i%6==0){
                PetscPrintf(PETSC_COMM_WORLD,"\n***                  ");
            }
            PetscPrintf(PETSC_COMM_WORLD,"%8d ",GetIthElmtJthConn(e,i));
        }
        if(GetIthElmtNodesNum(e)<5){
            for(PetscInt i=1;i<=2*(5-GetIthElmtNodesNum(e))-1;++i){
                if(GetIthElmtNodesNum(e)==4) PetscPrintf(PETSC_COMM_WORLD,"           ");
                if(GetIthElmtNodesNum(e)==3) PetscPrintf(PETSC_COMM_WORLD,"      ");
                if(GetIthElmtNodesNum(e)==2) PetscPrintf(PETSC_COMM_WORLD,"     ");
                if(GetIthElmtNodesNum(e)==1) PetscPrintf(PETSC_COMM_WORLD,"     ");
            }
            if(GetIthElmtNodesNum(e)==3) PetscPrintf(PETSC_COMM_WORLD,"  ");
            if(GetIthElmtNodesNum(e)==2) PetscPrintf(PETSC_COMM_WORLD,"    ");
            if(GetIthElmtNodesNum(e)==1) PetscPrintf(PETSC_COMM_WORLD,"   ");
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"  ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"***\n");
    }
    

    PetscPrintf(PETSC_COMM_WORLD,"***  node coornidates (node id, x, y, z and weight)                   ***\n");
    for(PetscInt i=1;i<=GetNodesNum();++i){
        PetscPrintf(PETSC_COMM_WORLD,"***  %9d:%13.4e,%13.4e,%13.4e, %11.3e ***\n",i,GetIthNodeJthCoord(i,1),
        GetIthNodeJthCoord(i,2),GetIthNodeJthCoord(i,3),GetIthNodeJthCoord(i,0));
    }
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
}