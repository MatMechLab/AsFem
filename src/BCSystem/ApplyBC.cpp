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
//+++ Purpose: here we apply the boundary conditions we defined in
//+++          each [bcs] sub blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"
#include "DofHandler/DofHandler.h"

void BCSystem::ApplyBC(const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const FECalcType &calctype,const double &t,const double (&ctan)[2],Vec &U,Mat &AMATRIX,Vec &RHS){
    double bcvalue;
    vector<string> bcnamelist;
    int DofIndex;
    if(ctan[0]){}
    for(auto it:_BCBlockList){
        bcvalue=it._BCValue;
        if(it._IsTimeDependent) bcvalue=it._BCValue*t;
        DofIndex=it._DofID;
        bcnamelist=it._BoundaryNameList;
        if(it._BCType==BCType::DIRICHLETBC){
            ApplyDirichletBC(mesh,dofHandler,calctype,DofIndex,bcvalue,bcnamelist,U,AMATRIX,RHS);
        }
        else if(it._BCType==BCType::NEUMANNBC){
            if(calctype==FECalcType::ComputeResidual){
                ApplyNeumannBC(mesh,dofHandler,fe,DofIndex,bcvalue,bcnamelist,RHS);
            }
        }
    }
}
//****************************************************
void BCSystem::ApplyInitialBC(const Mesh &mesh,const DofHandler &dofHandler,const double &t,Vec &U){
    PetscReal bcvalue;
    vector<string> bclist;
    string bcname;
    PetscInt DofIndex;
    PetscInt i,j,e,ee,iInd;
    int rankne,eStart,eEnd;

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    
    for(auto it:_BCBlockList){
        bcvalue=it._BCValue;
        if(it._IsTimeDependent) bcvalue=t*it._BCValue;
        bclist=it._BoundaryNameList;
        DofIndex=it._DofID;
        if(it._BCType==BCType::DIRICHLETBC){
            for(unsigned int ibc=0;ibc<bclist.size();++ibc){
                bcname=bclist[ibc];
                rankne=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname)/_size;
                eStart=_rank*rankne;
                eEnd=(_rank+1)*rankne;
                if(_rank==_size-1) eEnd=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname);
                for(e=eStart;e<eEnd;++e){
                    ee=mesh.GetBulkMeshIthElmtIDViaPhyName(bcname,e+1);
                    // cout<<"bcname="<<bcname<<", bcvalue="<<bcvalue<<", e="<<ee<<":"<<endl;
                    for(i=1;i<=mesh.GetBulkMeshIthElmtNodesNumViaPhyName(bcname,e+1);++i){
                        j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                        VecSetValues(U,1,&iInd,&bcvalue,INSERT_VALUES);
                        // cout<<iInd+1<<" ";
                    }
                    // cout<<endl;
                }
            }
        }
        else{
            continue;
        }
    }

    VecAssemblyBegin(U);
    VecAssemblyEnd(U);
}