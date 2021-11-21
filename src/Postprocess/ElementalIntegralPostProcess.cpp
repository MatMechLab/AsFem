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
//+++ Purpose: this pps can calculate the intergrated value over
//+++          specific domain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocess.h"

double Postprocess::ElementalIntegralPostProcess(vector<string> domainnamelist,string variablename,
                                                 const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const SolutionSystem &solutionSystem){
    double value=0.0,dofvalue;
    int nDim,nNodesPerElmt;
    int i,j,e,ee,iInd,gpInd,DofIndex;
    Nodes elNodes;
    elNodes.InitNodes(27);
    double xi,eta,zeta,JxW;
    double elU[27];


    DofIndex=dofHandler.GetDofIDviaDofName(variablename);

    if(domainnamelist.size()<1){
        MessagePrinter::PrintErrorTxt("error detected in ElementalIntegralPostProcess, we can not find any domain name, "
                                      "please check your input file or your mesh");
        MessagePrinter::AsFem_Exit();
    }
    if(DofIndex<1){
        MessagePrinter::PrintErrorTxt("error detected in ElementalIntegralPostProcess,"
                                      "we can not find DoF(name="+variablename+")"
                                      +", please check either your input file");
        MessagePrinter::AsFem_Exit();
    }

    VecScatterCreateToAll(solutionSystem._Unew,&_scatteru,&_Useq);
    VecScatterBegin(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    value=0.0;

    for(const auto &domainname:domainnamelist){
        for(e=1;e<=mesh.GetBulkMeshElmtsNumViaPhysicalName(domainname);e++){
            nDim=mesh.GetBulkMeshDimViaPhyName(domainname);
            if(nDim!=mesh.GetBulkMeshDim()){
                MessagePrinter::PrintErrorTxt("error detected in ElementalIntegralPostProcess,"
                                              "your mesh dim="+to_string(mesh.GetBulkMeshDim())+", however your element dim="
                                              +to_string(nDim)+", they are not match");
                MessagePrinter::AsFem_Exit();
            }
            nNodesPerElmt=mesh.GetBulkMeshNodesNumPerElmtViaPhysicalName(domainname);
            ee=mesh.GetBulkMeshIthElmtIDViaPhyName(domainname,e);
            // get the dof value for each nodal point
            for(i=1;i<=nNodesPerElmt;++i){
                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                iInd=dofHandler.GetBulkMeshIthNodeJthDofIndex(j,DofIndex)-1;
                VecGetValues(_Useq,1,&iInd,&dofvalue);
                elU[i-1]=dofvalue;
            }
            if(nDim==0){
                MessagePrinter::PrintErrorTxt("you can not get the 'elemental value' of a point, the dimension for elemental integration postprocess must be 1, 2 or 3");
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
                    // now we can do the gauss point integration
                    for(i=1;i<=nNodesPerElmt;++i){
                        value+=fe._BulkShp.shape_value(i)*elU[i-1]*JxW;
                    }
                }
            }
        }
    }
    // delete scatter
    VecScatterDestroy(&_scatteru);
    VecDestroy(&_Useq);

    return value;
}
