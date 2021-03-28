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
//+++ Date   : 2021.02.21
//+++ Purpose: here we apply nodal type dirichlet boundary condition
//             via the penalty method
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "BCSystem/BCSystem.h"
#include "DofHandler/DofHandler.h"

void BCSystem::ApplyNodalDirichletBC(const Mesh &mesh,const DofHandler &dofHandler,const FECalcType &calctype,const int &dofindex,const double &bcvalue,const vector<string> &bcnamelist,Vec &U,Mat &K,Vec &RHS){
    PetscInt nodeid,e;
    PetscInt iInd;
    const PetscScalar fix=0.0;
    int rankne,eStart,eEnd;

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    for(auto bcname:bcnamelist){
        rankne=mesh.GetBulkMeshNodeIDsNumViaPhysicalName(bcname)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetBulkMeshNodeIDsNumViaPhysicalName(bcname);
        for(e=eStart;e<eEnd;++e){
            nodeid=mesh.GetBulkMeshIthNodeIDViaPhyName(bcname,e+1);
            iInd=dofHandler.GetIthNodeJthDofIndex(nodeid,dofindex)-1;
            if(calctype==FECalcType::ComputeResidual) {
                VecSetValues(RHS,1,&iInd,&fix,INSERT_VALUES);
            }
            else if(calctype==FECalcType::ComputeJacobian){
                MatSetValues(K,1,&iInd,1,&iInd,&_PenaltyFactor,INSERT_VALUES);
            }
            VecSetValues(U,1,&iInd,&bcvalue,INSERT_VALUES);
        }// end-of-nodes-loop
    }// end-of-boundary-name-loop
}
