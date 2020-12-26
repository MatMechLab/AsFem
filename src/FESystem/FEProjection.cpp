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
//+++ Date   : 2020.12.25
//+++ Purpose: assemble and project gauss points' quantities to nodal
//+++          point
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/FESystem.h"

void FESystem::AssembleLocalProjectionToGlobal(const int &iInd,const int &nproj,
                                            const double &detjac,const double &test,
                                            const vector<double> &projvec,Vec &ProjVec){
    double w;
    int k,jInd;
    jInd=iInd*(nproj+1)+0;
    w=detjac*test;
    VecSetValue(ProjVec,jInd,w,ADD_VALUES);
    for(k=1;k<=nproj;k++){
        jInd=iInd*(nproj+1)+k;
        VecSetValue(ProjVec,jInd,w*projvec[k-1],ADD_VALUES);
    }
}
//******************************************************
    
//******************************************************
void FESystem::Projection(const int &nTotalNodes,const int &nproj,Vec &ProjVec){
    int i,j,iInd;
    double value,weight,newvalue;
    VecAssemblyBegin(ProjVec);
    VecAssemblyEnd(ProjVec);
    Vec ProjCopy;
    VecDuplicate(ProjVec,&ProjCopy);
    VecScatterCreateToAll(ProjVec,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,ProjVec,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,ProjVec,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);

    VecSet(ProjCopy,0.0);

    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);

    int rankne=nTotalNodes/_size;
    int eStart=_rank*rankne;
    int eEnd=(_rank+1)*rankne;
    if(_rank==_size-1) eEnd=nTotalNodes;

    for(int ee=eStart;ee<eEnd;++ee){
        i=ee+1;
        iInd=(i-1)*(1+nproj)+0;
        VecGetValues(_ProjSeq,1,&iInd,&weight);
        for(j=1;j<=nproj;j++){
            iInd=(i-1)*nproj+j;
            VecGetValues(_ProjSeq,1,&iInd,&value);
            VecSetValue(ProjCopy,iInd,value,ADD_VALUES);
            if(abs(value/weight)>1.0e-14){
                newvalue=value/weight;
            }
            else{
                newvalue=1.0e-14;
            }
            VecSetValue(ProjCopy,iInd,newvalue,ADD_VALUES);
        }
    }
    VecAssemblyBegin(ProjCopy);
    VecAssemblyEnd(ProjCopy);

    VecCopy(ProjCopy,ProjVec);
    
    VecScatterDestroy(&_scatterproj);
    VecDestroy(&_ProjSeq);
}