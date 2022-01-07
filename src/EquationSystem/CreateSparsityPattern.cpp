//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.26
//+++ Purpose: here we create the sparsity pattern of our system 
//+++          matrix according to our dof maps
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "EquationSystem/EquationSystem.h"

void EquationSystem::CreateSparsityPattern(DofHandler &dofHandler){
    PetscMPIInt rank,size;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    int rankne=dofHandler.GetBulkMeshBulkElmtsNum()/size;
    int eStart=rank*rankne;
    int eEnd=(rank+1)*rankne;
    if(rank==size-1) eEnd=dofHandler.GetBulkMeshBulkElmtsNum();
    int nDofs;
    vector<double> localK;
    vector<int> elDofs;

    nDofs=dofHandler.GetMaxDofsNumPerBulkElmt();
    elDofs.resize(nDofs,0);
    localK.resize(nDofs*nDofs,0.0);
    
    for(int e=eStart;e<eEnd;++e){
        dofHandler.GetBulkMeshIthBulkElmtDofIndex0(e+1,elDofs,localK);
        nDofs=dofHandler.GetBulkMeshIthBulkElmtDofsNum(e+1);
        MatSetValues(_AMATRIX,nDofs,elDofs.data(),nDofs,elDofs.data(),localK.data(),ADD_VALUES);
    }
    MatAssemblyBegin(_AMATRIX,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_AMATRIX,MAT_FINAL_ASSEMBLY);
    localK.resize(0);
    elDofs.resize(0);
    //*****************************************************************************
    //*** once the sparsity pattern is ready, we need to disable the new allocation 
    //*** to our matrix
    //*****************************************************************************
    MatSetOption(_AMATRIX,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
}
