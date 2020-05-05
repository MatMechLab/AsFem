//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "BCs/BCSystem.h"

void BCSystem::ApplyBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const PetscReal &t,const PetscReal (&ctan)[2],Mat &K,Vec &RHS,Vec &U){
    //*********************************
    //*** just to get rid of warning
    //*********************************
    if(fe.GetDim()){};
    if(ctan[0]){};
    
    PetscReal bcvalue;
    vector<string> bclist;
    PetscInt DofIndex;
    
    for(auto it:_BCBlockList){
        bcvalue=it._BCValue;
        if(it._IsTimeDependent) bcvalue=t*it._BCValue;
        bclist=it._BoundaryNameList;
        DofIndex=it._DofIndex;
        if(it._BCType==BCType::DirichletBC){
            ApplyDirichletBC(mesh,dofHandler,DofIndex,bcvalue,bclist,K,RHS,U);
        }
        else if(it._BCType==BCType::NeumannBC){
            ApplyNeumannBC(mesh,dofHandler,fe,DofIndex,bcvalue,bclist,RHS);
        }
        else if(it._BCType==BCType::PresetBC){
            ApplyPresetBC(mesh,dofHandler,fe,DofIndex,bcvalue,bclist,ctan,K,RHS,U);
        }
        else if(it._BCType==BCType::PressureBC){
            // cout<<"t="<<t<<", bcvalue="<<bcvalue<<endl;
            ApplyPressureBC(mesh,dofHandler,fe,DofIndex,bcvalue,bclist,RHS);
        }
        else if(it._BCType==BCType::User1BC){
            ApplyUser1BC(mesh,dofHandler,fe,DofIndex,bcvalue,bclist,ctan,K,RHS,U);
        }
        else{
            continue;
        }
    }
    VecAssemblyBegin(RHS);
    VecAssemblyEnd(RHS);
    VecAssemblyBegin(U);
    VecAssemblyEnd(U);
    MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
}
//***************************************
void BCSystem::ApplyInitialBC(Mesh &mesh,DofHandler &dofHandler,const PetscReal &t,Vec &U){
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