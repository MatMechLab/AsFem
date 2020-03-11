//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FESystem/FESystem.h"

void FESystem::Projection(const int &nNodes,Vec &Proj,const vector<double> &elProj,ShapeFun &shp,const double &DetJac){
    double w;
    int j,k,iInd,jInd;
    for(j=1;j<=nNodes;++j){
        iInd=_elConn[j-1]-1;
        jInd=iInd*(_nProj+1)+0;
        w=DetJac*shp.shape_value(j);
        VecSetValue(Proj,jInd,w,ADD_VALUES);
        for(k=1;k<=_nProj;++k){
            jInd=iInd*(_nProj+1)+k;
            VecSetValue(Proj,jInd,w*elProj[k-1],ADD_VALUES);
        }
    }
}