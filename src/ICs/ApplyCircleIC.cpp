//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ICs/ICSystem.h"

void ICSystem::ApplyCircleIC(const vector<double> Params,const int &DofIndex,Mesh &mesh,DofHandler &dofHandler,Vec &U){
    if(Params.size()<6){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for circle IC, you need at least 5 params            !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        x0,y0,radius,rwidth,valuein,valueout are expected    !!!   ***\n");
        Msg_AsFem_Exit();
    }


    int iInd,i,ii;
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    int rankn=mesh.GetNodesNum()/_size;
    int iStart=_rank*rankn;
    int iEnd=(_rank+1)*rankn;
    if(_rank==_size-1) iEnd=mesh.GetNodesNum();
    PetscScalar value;

    double x,y,x0,y0,radius,rwidth,theta,dist,val1,val2;
    const double PI=3.14159265359;

    x0=Params[0];y0=Params[1];
    radius=Params[2];rwidth=Params[3];
    val1=Params[4];val2=Params[5];

    for(ii=iStart;ii<iEnd;++ii){
        i=ii+1;
        x=mesh.GetIthNodeJthCoord(i,1);
        y=mesh.GetIthNodeJthCoord(i,2);
        dist=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
        if(dist<=radius-rwidth/2.0){
            value=val1;
        }
        else if(dist<radius+rwidth/2.0 && dist>radius-rwidth/2.0){
            theta=(dist-radius+rwidth/2.0)/rwidth;
            value=val2+(val1-val2)*(1+cos(theta*PI))/2.0;
        }
        else{
            value=val2;
        }
        iInd=dofHandler.GetIthNodeJthDofIndex(i,DofIndex)-1;
        VecSetValue(U,iInd,value,INSERT_VALUES);
    }
}