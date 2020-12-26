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

void BCSystem::ApplyBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const FECalcType &calctype,const double &t,const double (&ctan)[2],Vec &U,Mat &AMATRIX,Vec &RHS){
    double bcvalue;
    vector<string> bcnamelist;
    int DofIndex;
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
                ApplyNeumannBC();
            }
        }
    }
}
//****************************************************
void BCSystem::ApplyInitialBC(Mesh &mesh,DofHandler &dofHandler,const double &t,Vec &U){
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
        DofIndex=it._DofIndex;
        if(it._BCType==BCType::DirichletBC){
            for(unsigned int ibc=0;ibc<bclist.size();++ibc){
                bcname=bclist[ibc];
                rankne=mesh.GetElmtsNumViaPhyName(bcname)/_size;
                eStart=_rank*rankne;
                eEnd=(_rank+1)*rankne;
                if(_rank==_size-1) eEnd=mesh.GetElmtsNumViaPhyName(bcname);
                for(e=eStart;e<eEnd;++e){
                    ee=mesh.GetIthElmtIndexViaPhyName(bcname,e+1);
                    // cout<<"bcname="<<bcname<<", bcvalue="<<bcvalue<<", e="<<ee<<":"<<endl;
                    for(i=1;i<=mesh.GetIthElmtNodesNumViaPhyName(bcname,e+1);++i){
                        j=mesh.GetIthElmtJthConn(ee,i);
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