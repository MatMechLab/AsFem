//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "EquationSystem/EquationSystem.h"

void EquationSystem::CreateSparsityPattern(DofHandler &dofHandler){
    PetscMPIInt rank,size;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    // PetscPrintf(PETSC_COMM_WORLD,"rank=%2d:",rank);

    int rankne=dofHandler.GetBulkElmtsNum()/size;
    int eStart=rank*rankne;
    int eEnd=(rank+1)*rankne;
    if(rank==size-1) eEnd=dofHandler.GetBulkElmtsNum();
    int nDofs;
    PetscReal localK[270];
    int conn[27];
    for(int i=0;i<270;++i){localK[i]=0.0;}

    for(int e=eStart;e<eEnd;++e){
        dofHandler.GetIthBulkElmtDofIndex0(e+1,conn);
        nDofs=dofHandler.GetIthBulkElmtDofsNum(e+1);
        MatSetValues(_AMATRIX,nDofs,conn,nDofs,conn,localK,ADD_VALUES);
        // for(int i=0;i<nDofs;++i){
        //     PetscPrintf(PETSC_COMM_WORLD,"%5d ",conn[i]);
        // }
        // PetscPrintf(PETSC_COMM_WORLD,"\n");
    }
    MatAssemblyBegin(_AMATRIX,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_AMATRIX,MAT_FINAL_ASSEMBLY);

    MatSetOption(_AMATRIX,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
    
}