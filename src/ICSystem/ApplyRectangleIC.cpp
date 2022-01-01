//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.30
//+++ Purpose: Apply the cubic IC to U vector
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

void ICSystem::ApplyRectangleIC(const int &DofIndex, const vector<double> &Parameters, const vector<string> &DomainList,
                                const Mesh &mesh, const DofHandler &dofHandler, Vec &U) {
    if(Parameters.size()<5){
        MessagePrinter::PrintErrorTxt("for constant IC, you need at least four parameter, params='x0 y0 dx dy value' is expected in [ics]");
        MessagePrinter::AsFem_Exit();
    }

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    double x0,y0,dx,dy,value;
    double x,y;
    x0=Parameters[0];y0=Parameters[1];
    dx=Parameters[2];dy=Parameters[3];
    value=Parameters[6];

    int rankne,eStart,eEnd,e,ee,i,j,iInd;
    for(auto domain:DomainList){
        rankne=mesh.GetBulkMeshElmtsNumViaPhysicalName(domain)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetBulkMeshElmtsNumViaPhysicalName(domain);
        for(e=eStart;e<eEnd;++e){
            ee=mesh.GetBulkMeshIthElmtIDViaPhyName(domain,e+1);//global id
            for(i=1;i<=mesh.GetBulkMeshIthElmtNodesNum(ee);++i){
                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                x=mesh.GetBulkMeshIthNodeJthCoord(j,1);
                y=mesh.GetBulkMeshIthNodeJthCoord(j,2);
                if((x>=x0-dx && x<=x0+dx)&&(y>=y0-dy&& y<=y0+dy)){
                    iInd=dofHandler.GetBulkMeshIthNodeJthDofIndex(j,DofIndex)-1;
                    VecSetValues(U,1,&iInd,&value,INSERT_VALUES);
                }
            }
        }
    }
}
