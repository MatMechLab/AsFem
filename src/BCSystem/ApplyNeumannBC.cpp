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
//+++ Purpose: here we apply the neumann boundary condition to our
//+++          residual
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"
#include "DofHandler/DofHandler.h"

void BCSystem::ApplyNeumannBC(const Mesh &mesh,const DofHandler &dofHandler,FE &fe,
                        const int &DofIndex,const double &bcvalue,const vector<string> &bcnamelist,
                        Vec &RHS){
    PetscInt i,j,e,ee,gpInd;
    PetscInt iInd;
    PetscScalar value;
    string bcname;
    int rankne,eStart,eEnd;


    for(auto bcname:bcnamelist){
        MPI_Comm_size(PETSC_COMM_WORLD,&_size);
        MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

        rankne=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname);

        _nDim=mesh.GetBulkMeshDimViaPhyName(bcname);
        _nNodesPerBCElmt=mesh.GetBulkMeshNodesNumPerElmtViaPhysicalName(bcname);

        for(e=eStart;e<eEnd;++e){
            ee=mesh.GetBulkMeshIthElmtIDViaPhyName(bcname,e+1);
            if(_nDim==0){
                // for point case,(bulk dim=1, bc dim=0)
                for(i=1;i<=_nNodesPerBCElmt;++i){
                    j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                    VecSetValue(RHS,iInd,bcvalue,ADD_VALUES);
                }
            }
            else if(_nDim==1){
                // for line case (bulk dim>=2, bc dim=1)
                mesh.GetBulkMeshIthElmtNodes(ee,_elNodes);
                for(gpInd=1;gpInd<=fe._LineQPoint.GetQpPointsNum();++gpInd){
                    _xi=fe._LineQPoint(gpInd,1);
                    fe._LineShp.Calc(_xi,_elNodes,true);
                    _JxW=fe._LineShp.GetDetJac()*fe._LineQPoint(gpInd,0);
                    for(i=1;i<=_nNodesPerBCElmt;++i){
                        j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                        value=fe._LineShp.shape_value(i)*bcvalue*_JxW;
                        VecSetValue(RHS,iInd,value,ADD_VALUES);
                    }
                }
            }
            else if(_nDim==2){
                // for surface case (bulk dim=3, bc dim=2)
                mesh.GetBulkMeshIthElmtNodes(ee,_elNodes);
                for(gpInd=1;gpInd<=fe._SurfaceQPoint.GetQpPointsNum();++gpInd){
                    _xi=fe._SurfaceQPoint(gpInd,1);
                    _eta=fe._SurfaceQPoint(gpInd,2);
                    fe._SurfaceShp.Calc(_xi,_eta,_elNodes,true);
                    _JxW=fe._SurfaceShp.GetDetJac()*fe._SurfaceQPoint(gpInd,0);
                    for(i=1;i<=_nNodesPerBCElmt;++i){
                        j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                        value=fe._SurfaceShp.shape_value(i)*bcvalue*_JxW;
                        VecSetValue(RHS,iInd,value,ADD_VALUES);
                    }
                }
            }
        }
    }
}