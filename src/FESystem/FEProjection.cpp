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
    VecScatterCreateToAll(ProjVec,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,ProjVec,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,ProjVec,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    for(i=1;i<=nTotalNodes;i++){
        iInd=(i-1)*(1+nproj)+0;
        VecGetValues(_ProjSeq,1,&iInd,&weight);
        for(j=1;j<=nproj;j++){
            iInd=(i-1)*nproj+j;
            VecGetValues(_ProjSeq,1,&iInd,&value);
            if(abs(value/weight)>1.0e-14){
                newvalue=value/weight;
            }
            else{
                newvalue=1.0e-14;
            }
            VecSetValue(ProjVec,iInd,newvalue,INSERT_VALUES);
        }
    }
    VecAssemblyBegin(ProjVec);
    VecAssemblyEnd(ProjVec);
    
    VecScatterDestroy(&_scatterproj);
    VecDestroy(&_ProjSeq);
}