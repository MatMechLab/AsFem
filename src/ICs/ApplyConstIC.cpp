//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ICs/ICSystem.h"

void ICSystem::ApplyConstIC(const vector<double> Params,const int &DofIndex,Mesh &mesh,DofHandler &dofHandler,Vec &U){
    if(Params.size()<1){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for constic, you need at least one param             !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        params=constvalue is expected in [ics]               !!!   ***\n");
        Msg_AsFem_Exit();
    }
    int iInd,i,ii;
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    int rankn=mesh.GetNodesNum()/_size;
    int iStart=_rank*rankn;
    int iEnd=(_rank+1)*rankn;
    if(_rank==_size-1) iEnd=mesh.GetNodesNum();

    for(ii=iStart;ii<iEnd;++ii){
        i=ii+1;
        iInd=dofHandler.GetIthNodeJthDofIndex(i,DofIndex)-1;
        VecSetValue(U,iInd,Params[0],INSERT_VALUES);
    }
}