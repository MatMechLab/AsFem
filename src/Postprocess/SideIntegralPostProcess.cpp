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
//+++ Date   : 2021.03.22
//+++ Purpose: this pps can calculate the intergrated value over
//+++          specific side-set for dof
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocess.h"

double Postprocess::SideIntegralPostProcess(vector<string> sidenamelist,string dofname,
                                            const Mesh &mesh,const DofHandler &dofHandler,FE &fe,
                                            const SolutionSystem &solutionSystem){
    double dofvalue,value;
    int nDim,nNodesPerElmt;
    int i,j,e,ee,iInd,gpInd,DofIndex;
    Nodes elNodes;
    elNodes.InitNodes(27);
    double xi,eta,JxW;
    double elU[27];

    VecScatterCreateToAll(solutionSystem._Unew,&_scatteru,&_Useq);
    VecScatterBegin(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    DofIndex=dofHandler.GetDofIDviaDofName(dofname);
    value=0.0;

    if(DofIndex<1){
        MessagePrinter::PrintErrorTxt("error detected in SideIntegralPostProcess,"
                                      "we can not find dof(name="+dofname+")"
                                      +", please check either your input file");
        MessagePrinter::AsFem_Exit();
    }
    if(sidenamelist.size()<1){
        MessagePrinter::PrintErrorTxt("error detected in SideIntegralPostProcess,"
                                      "we can not find any side sets , please check either your input file");
        MessagePrinter::AsFem_Exit();
    }
    for(const auto &sidename:sidenamelist){
        for(e=1;e<=mesh.GetBulkMeshElmtsNumViaPhysicalName(sidename);e++){
            nDim=mesh.GetBulkMeshDimViaPhyName(sidename);
            if(nDim==mesh.GetBulkMeshDim()){
                MessagePrinter::PrintErrorTxt("error detected in ProjVariableSideIntegralPostProcess,"
                                              "your mesh dim="+to_string(mesh.GetBulkMeshDim())+", however your side dim="
                                              +to_string(nDim)+" is equal to your bulk mesh");
                MessagePrinter::AsFem_Exit();
            }
            nNodesPerElmt=mesh.GetBulkMeshNodesNumPerElmtViaPhysicalName(sidename);
            ee=mesh.GetBulkMeshIthElmtIDViaPhyName(sidename,e);
            // get the dof value for each nodal point
            for(i=1;i<=nNodesPerElmt;++i){
                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                VecGetValues(_Useq,1,&iInd,&dofvalue);
                elU[i-1]=dofvalue;
            }
            if(nDim==0){
                // the bulk mesh is 1d line, so the side-set is nodal point(dim=0)
                value+=dofvalue;
            }
            else if(nDim==1){
                // for the 1d-line(dim=1) in 2d domain(dim=2)
                mesh.GetBulkMeshIthElmtNodes(ee,elNodes);
                for(gpInd=1;gpInd<=fe._LineQPoint.GetQpPointsNum();++gpInd){
                    xi=fe._LineQPoint(gpInd,1);
                    fe._LineShp.Calc(xi,elNodes,true);
                    JxW=fe._LineShp.GetDetJac()*fe._LineQPoint(gpInd,0);
                    // now we can do the gauss point integration
                    for(i=1;i<=nNodesPerElmt;++i){
                        value+=fe._LineShp.shape_value(i)*elU[i-1]*JxW;
                    }
                }
            }
            else if(nDim==2){
                // for the 2d-surface(dim=2) in 3d domain(dim=3)
                mesh.GetBulkMeshIthElmtNodes(ee,elNodes);
                for(gpInd=1;gpInd<=fe._SurfaceQPoint.GetQpPointsNum();++gpInd){
                    xi=fe._SurfaceQPoint(gpInd,1);
                    eta=fe._SurfaceQPoint(gpInd,2);
                    fe._SurfaceShp.Calc(xi,eta,elNodes,true);
                    JxW=fe._SurfaceShp.GetDetJac()*fe._SurfaceQPoint(gpInd,0);
                    // now we can do the gauss point integration
                    for(i=1;i<=nNodesPerElmt;++i){
                        value+=fe._SurfaceShp.shape_value(i)*elU[i-1]*JxW;
                    }
                }
            }
            else{
                MessagePrinter::PrintErrorTxt("error detected in ProjVariableSideIntegralPostProcess,"
                                              "your mesh dim="+to_string(mesh.GetBulkMeshDim())+", however your side dim="
                                              +to_string(nDim)+" is equal to your bulk mesh");
                MessagePrinter::AsFem_Exit();
            }
        }
    }
    // delete scatter
    VecScatterDestroy(&_scatteru);
    VecDestroy(&_Useq);

    return value;
}