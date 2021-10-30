//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.30
//+++ Purpose: Apply the random IC to U vector
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

void ICSystem::ApplyRandomIC(const int &DofIndex, const vector<double> &Parameters, const vector<string> &DomainList,
                             const Mesh &mesh, const DofHandler &dofHandler, Vec &U){
    if(Parameters.size()<2){
        MessagePrinter::PrintErrorTxt("for random IC, you need at least two parameters, params='low-value high-value' is expected in [ics]");
        MessagePrinter::AsFem_Exit();
    }

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    PetscRandomCreate(PETSC_COMM_WORLD,&_rnd);
    PetscRandomSetInterval(_rnd,Parameters[0],Parameters[1]);
    PetscRandomSetType(_rnd,PETSCRAND);

    int rankne,eStart,eEnd,e,ee,i,j,iInd;
    double value;
    for(auto domain:DomainList){
        rankne=mesh.GetBulkMeshElmtsNumViaPhysicalName(domain)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetBulkMeshElmtsNumViaPhysicalName(domain);
        for(e=eStart;e<eEnd;++e){
            ee=mesh.GetBulkMeshIthElmtIDViaPhyName(domain,e+1);//global id
            for(i=1;i<=mesh.GetBulkMeshIthElmtNodesNum(ee);++i){
                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                PetscRandomGetValue(_rnd,&value);
                VecSetValues(U,1,&iInd,&value,INSERT_VALUES);
            }
        }
    }

    PetscRandomDestroy(&_rnd);
}
