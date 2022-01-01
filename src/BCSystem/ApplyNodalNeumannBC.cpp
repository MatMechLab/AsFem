//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.02.21
//+++ Purpose: here we apply the neumann boundary condition to
//+++          nodal sets
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"
#include "DofHandler/DofHandler.h"

void BCSystem::ApplyNodalNeumannBC(const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const vector<int> &dofsindex,const double &bcvalue,const vector<string> &bcnamelist,Vec &RHS){
    PetscInt j,e;
    PetscInt iInd;
    int rankne,eStart,eEnd;
    if(fe.GetDim()){}

    for(const auto &bcname:bcnamelist){
        MPI_Comm_size(PETSC_COMM_WORLD,&_size);
        MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
        rankne=mesh.GetBulkMeshNodeIDsNumViaPhysicalName(bcname)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetBulkMeshNodeIDsNumViaPhysicalName(bcname);
        for(e=eStart;e<eEnd;++e){
            j=mesh.GetBulkMeshIthNodeIDViaPhyName(bcname,e+1);
            for(const auto &id:dofsindex){
                iInd=dofHandler.GetBulkMeshIthNodeJthDofIndex(j,id)-1;
                VecSetValue(RHS,iInd,bcvalue,ADD_VALUES);
            }
        }
    }

}
