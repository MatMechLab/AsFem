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
//+++ Date   : 2020.12.26
//+++ Purpose: here we apply dirichlet boundary condition via the 
//+++          penalty method
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "BCSystem/BCSystem.h"
#include "DofHandler/DofHandler.h"

void BCSystem::ApplyDirichletBC(const FECalcType &calctype,const BCType &bctype,const vector<string> bcnamelist,const vector<int> &dofindex,const double &bcvalue,const vector<double> &params,const Mesh &mesh,const DofHandler &dofHandler,Vec &U,Mat &K,Vec &RHS){
    PetscInt i,j,k,e,ee;
    PetscInt iInd;
    int rankne,eStart,eEnd;
    vector<int> dofids;

    dofids.resize(dofindex.size()+1,0);
    _elmtinfo.nDofs=static_cast<int>(dofindex.size());


    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    for(auto bcname:bcnamelist){
        rankne=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname);

        _elmtinfo.nDim=mesh.GetBulkMeshDimViaPhyName(bcname);

        for(e=eStart;e<eEnd;++e){
            ee=mesh.GetBulkMeshIthElmtIDViaPhyName(bcname,e+1);//global id
            _elmtinfo.nNodes=mesh.GetBulkMeshIthElmtNodesNum(ee);
            for(i=1;i<=mesh.GetBulkMeshIthElmtNodesNum(ee);++i){
                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                _elmtinfo.gpCoords(1)=mesh.GetBulkMeshIthNodeJthCoord(j,1);
                _elmtinfo.gpCoords(2)=mesh.GetBulkMeshIthNodeJthCoord(j,2);
                _elmtinfo.gpCoords(3)=mesh.GetBulkMeshIthNodeJthCoord(j,3);
                for(k=1;k<=static_cast<int>(dofindex.size());k++){
                    iInd=dofHandler.GetIthNodeJthDofIndex(j,dofindex[k-1])-1;
                    dofids[k-1]=iInd;
                }

                switch (bctype) {
                    case BCType::DIRICHLETBC:
                        DirichletBC::ComputeBCValue(calctype,bcvalue,params,_elmtinfo,dofids,_elmtinfo.gpCoords,K,RHS,U);
                        break;
                    case BCType::USER1DIRICHLETBC:
                        User1DirichletBC::ComputeBCValue(calctype,bcvalue,params,_elmtinfo,dofids,_elmtinfo.gpCoords,K,RHS,U);
                        break; 
                    case BCType::USER2DIRICHLETBC:
                        User2DirichletBC::ComputeBCValue(calctype,bcvalue,params,_elmtinfo,dofids,_elmtinfo.gpCoords,K,RHS,U);
                        break; 
                    case BCType::USER3DIRICHLETBC:
                        User3DirichletBC::ComputeBCValue(calctype,bcvalue,params,_elmtinfo,dofids,_elmtinfo.gpCoords,K,RHS,U);
                        break; 
                    case BCType::USER4DIRICHLETBC:
                        User4DirichletBC::ComputeBCValue(calctype,bcvalue,params,_elmtinfo,dofids,_elmtinfo.gpCoords,K,RHS,U);
                        break; 
                    case BCType::USER5DIRICHLETBC:
                        User5DirichletBC::ComputeBCValue(calctype,bcvalue,params,_elmtinfo,dofids,_elmtinfo.gpCoords,K,RHS,U);
                        break; 
                    default:
                        MessagePrinter::PrintErrorTxt("unsupported boundary condition type in ApplyDirichletBC, please check your code");
                        MessagePrinter::AsFem_Exit();
                        break;
                }
            }
        }
    }
}
