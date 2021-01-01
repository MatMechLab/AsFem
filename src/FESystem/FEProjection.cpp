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
//+++ Date   : 2020.12.25
//+++ Purpose: assemble and project gauss points' quantities to nodal
//+++          point
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/FESystem.h"

void FESystem::AssembleLocalProjectionToGlobal(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                            const vector<double> &elProj,Vec &ProjVec){
    double w;
    int j,k,jInd,iInd;
    for(j=1;j<=nNodes;j++){
        iInd=_elConn[j-1]-1;
        jInd=iInd*(_nProj+1)+0;
        w=DetJac*shp.shape_value(j);
        VecSetValue(ProjVec,jInd,w,ADD_VALUES);
        for(k=1;k<=_nProj;k++){
            jInd=iInd*(_nProj+1)+k;
            VecSetValue(ProjVec,jInd,w*elProj[k-1],ADD_VALUES);
        }
    }
}
//******************************************************
    
//******************************************************
void FESystem::Projection(const int &nTotalNodes,const int &nproj,Vec &ProjVec){
    int i,j,iInd;
    double value,weight,newvalue;

    VecAssemblyBegin(ProjVec);
    VecAssemblyEnd(ProjVec);

    VecScatterCreateToAll(ProjVec,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,ProjVec,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,ProjVec,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);

    VecSet(ProjVec,0.0);

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
        VecSetValue(ProjVec,iInd,weight,ADD_VALUES);
        for(j=1;j<=nproj;j++){
            iInd=(i-1)*(nproj+1)+j;
            VecGetValues(_ProjSeq,1,&iInd,&value);
            if(abs(weight)>1.0e-15){
                newvalue=value/weight;
            }
            else{
                newvalue=value;
            }
            VecSetValue(ProjVec,iInd,newvalue,ADD_VALUES);
        }
    }
    VecAssemblyBegin(ProjVec);
    VecAssemblyEnd(ProjVec);

    
    VecScatterDestroy(&_scatterproj);
    VecDestroy(&_ProjSeq);
}