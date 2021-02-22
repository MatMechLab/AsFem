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
//+++ Date   : 2021.02.22
//+++ Purpose: this pps can calculate the area of specific side or
//+++          boundary
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocess.h"

double Postprocess::VolumePostProcess(vector<string> domainnamelist,const Mesh &mesh,FE &fe){
    double volume=0.0;
    int nDim,nNodesPerElmt;
    int i,e,ee,gpInd;
    Nodes elNodes;
    elNodes.InitNodes(27);
    double xi,eta,zeta,JxW;


    volume=0.0;
    for(const auto &domainname:domainnamelist){
        for(e=1;e<=mesh.GetBulkMeshElmtsNumViaPhysicalName(domainname);e++){
            nDim=mesh.GetBulkMeshDimViaPhyName(domainname);
            if(nDim!=mesh.GetBulkMeshDim()){
                MessagePrinter::PrintErrorTxt("error detected in VolumePostProcess,"
                                              "your mesh dim="+to_string(mesh.GetBulkMeshDim())+", however your element dim="
                                              +to_string(nDim)+", they are not match");
                MessagePrinter::AsFem_Exit();
            }
            nNodesPerElmt=mesh.GetBulkMeshNodesNumPerElmtViaPhysicalName(domainname);
            ee=mesh.GetBulkMeshIthElmtIDViaPhyName(domainname,e);
            if(nDim==0){
                MessagePrinter::PrintErrorTxt("you can not get the 'volume' of a point, the dimension for volume postprocess must be 1, 2 or 3");
                MessagePrinter::AsFem_Exit();
            }
            else{
                // for 1D line
                mesh.GetBulkMeshIthElmtNodes(ee,elNodes);
                for(gpInd=1;gpInd<=fe._BulkQPoint.GetQpPointsNum();++gpInd){
                    xi=fe._BulkQPoint(gpInd,1);
                    if(nDim==1){
                        fe._BulkShp.Calc(xi,elNodes,true);
                    }
                    else if(nDim==2){
                        eta=fe._BulkQPoint(gpInd,2);
                        fe._BulkShp.Calc(xi,eta,elNodes,true);
                    }
                    else if(nDim==3){
                        eta=fe._BulkQPoint(gpInd,2);
                        zeta=fe._BulkQPoint(gpInd,3);
                        fe._BulkShp.Calc(xi,eta,zeta,elNodes,true);
                    }
                    JxW=fe._BulkShp.GetDetJac()*fe._BulkQPoint(gpInd,0);
                    for(i=1;i<=nNodesPerElmt;++i){
                        volume+=fe._BulkShp.shape_value(i)*1.0*JxW;
                    }
                }
            }
        }
    }
    return volume;
}