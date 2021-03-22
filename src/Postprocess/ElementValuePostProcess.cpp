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
//+++ Purpose: this pps can calculate the dofs value of specific elmtid
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocess.h"

double Postprocess::ElementValuePostProcess(const int &elmtid,string variablename,
                                         const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem){
    if(elmtid<1||elmtid>mesh.GetBulkMeshBulkElmtsNum()){
        MessagePrinter::PrintErrorTxt("elmtid="+to_string(elmtid)+" is invalid for ElementValuePostProcess");
        MessagePrinter::AsFem_Exit();
    }
    double elmtvalue=0.0,val;
    int DofIndex,j,iInd,i;
    DofIndex=dofHandler.GetDofIDviaDofName(variablename);
    if(DofIndex<1){
        MessagePrinter::PrintErrorTxt("error detected in ElementValuePostProcess,"
                                      "we can not find DoF(name="+variablename+")"
                                      +", please check either your input file");
        MessagePrinter::AsFem_Exit();
    }

    elmtvalue=0.0;

    VecScatterCreateToAll(solutionSystem._Unew,&_scatteru,&_Useq);
    VecScatterBegin(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    for(i=1;i<=mesh.GetBulkMeshIthBulkElmtNodesNum(elmtid);i++){
        j=mesh.GetBulkMeshIthBulkElmtJthNodeID(elmtid,i);
        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
        VecGetValues(_Useq,1,&iInd,&val);
        elmtvalue+=val;
    }
    // delete scatter
    VecScatterDestroy(&_scatteru);
    VecDestroy(&_Useq);

    return static_cast<double>(elmtvalue/mesh.GetBulkMeshIthBulkElmtNodesNum(elmtid));
}