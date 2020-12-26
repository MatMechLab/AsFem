//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.26
//+++ Purpose: here we apply dirichlet boundary condition via the 
//+++          penalty method
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "BCSystem/BCSystem.h"

void BCSystem::ApplyDirichletBC(Mesh &mesh,DofHandler &dofHandler,const FECalcType &calctype,const int &dofindex,const double &bcvalue,const vector<string> &bcnamelist,Vec &U,Mat &K,Vec &RHS){
    PetscInt i,j,e,ee;
    PetscInt iInd;
    const PetscScalar fix=0.0;
    int rankne,eStart,eEnd;

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    for(auto bcname:bcnamelist){
        rankne=mesh.GetElmtsNumViaPhysicalName(bcname)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetElmtsNumViaPhysicalName(bcname);

        for(e=eStart;e<eEnd;++e){
            ee=mesh.GetBulkMeshIthElmtIDViaPhyName(bcname,e+1);
            // cout<<"bcname="<<_bcname<<", bcvalue="<<_bcvalue<<", e="<<ee<<":";
            for(i=1;i<=mesh.GetBulkMeshIthElmtNodesNumViaPhyName(bcname,ee);++i){
                j=mesh.GetIthElmtJthNodeID(ee,i);
                iInd=dofHandler.GetIthNodeJthDofIndex(j,dofindex)-1;
                if(calctype==FECalcType::ComputeResidual) {
                    VecSetValues(RHS,1,&iInd,&fix,INSERT_VALUES);
                }
                else if(calctype==FECalcType::ComputeJacobian){
                    MatSetValues(K,1,&iInd,1,&iInd,&_PenaltyFactor,INSERT_VALUES);
                }
                VecSetValues(U,1,&iInd,&bcvalue,INSERT_VALUES);
            }
        }
    }
}
