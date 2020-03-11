//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "BCs/BCSystem.h"

void BCSystem::ApplyDirichletBC(Mesh &mesh,DofHandler &dofHandler,const PetscInt &DofIndex,const PetscReal &bcvalue,vector<string> bclist,Mat &K,Vec &RHS,Vec &U){
    PetscInt i,j,e,ee;
    PetscInt iInd;
    string bcname;
    const PetscScalar fix=0.0;
    int rankne,eStart,eEnd;

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    for(unsigned int ibc=0;ibc<bclist.size();++ibc){
        bcname=bclist[ibc];

        rankne=mesh.GetElmtsNumViaPhyName(bcname)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetElmtsNumViaPhyName(bcname);

        for(e=eStart;e<eEnd;++e){
            ee=mesh.GetIthElmtIndexViaPhyName(bcname,e+1);
            // cout<<"bcname="<<_bcname<<", bcvalue="<<_bcvalue<<", e="<<ee<<":";
            for(i=1;i<=mesh.GetIthElmtNodesNumViaPhyName(bcname,e+1);++i){
                j=mesh.GetIthElmtJthConn(ee,i);
                iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                VecSetValues(U,1,&iInd,&bcvalue,INSERT_VALUES);
                VecSetValues(RHS,1,&iInd,&fix,INSERT_VALUES);
                MatSetValues(K,1,&iInd,1,&iInd,&_PenaltyFactor,INSERT_VALUES);
                // cout<<iInd+1<<" ";
            }
            // cout<<endl;
        }
    }
}