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

void BCSystem::ApplyNodalDirichletBC(const FECalcType &calctype,const BCType &bctype,const vector<string> bcnamelist,const vector<int> &dofsindex,const double &bcvalue,const vector<double> &params,const Mesh &mesh,const DofHandler &dofHandler,Vec &U,Mat &K,Vec &RHS){
    PetscInt nodeid,e;
    PetscInt iInd;
    int rankne,eStart,eEnd;
    vector<int> dofsid;

    dofsid.resize(dofsindex.size(),0);
    _elmtinfo.nDofs=static_cast<int>(dofsindex.size());
    _elmtinfo.nDim=0;
    _elmtinfo.nNodes=1;
    
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    for(auto bcname:bcnamelist){
        rankne=mesh.GetBulkMeshNodeIDsNumViaPhysicalName(bcname)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetBulkMeshNodeIDsNumViaPhysicalName(bcname);
        for(e=eStart;e<eEnd;++e){
            nodeid=mesh.GetBulkMeshIthNodeIDViaPhyName(bcname,e+1);
            for(int i=0;i<_elmtinfo.nDofs;i++){
                iInd=dofHandler.GetBulkMeshIthNodeJthDofIndex(nodeid,dofsindex[i])-1;
                dofsid[i]=iInd;
            }
            switch (bctype) {
                case BCType::NODALDIRICHLETBC:
                    DirichletBC::ComputeBCValue(calctype,bcvalue,params,_elmtinfo,dofsid,_elmtinfo.gpCoords,K,RHS,U);
                    break;
                default:
                    MessagePrinter::PrintErrorTxt("unsupported boundary condition type in ApplyNodalDirichletBC, please check your code");
                    MessagePrinter::AsFem_Exit();
                    break;
                }
        }// end-of-nodes-loop
    }// end-of-boundary-name-loop
}
