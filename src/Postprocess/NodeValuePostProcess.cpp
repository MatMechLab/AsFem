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
//+++ Purpose: this pps can calculate the dofs value of specific nodeid
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocess.h"

double Postprocess::NodeValuePostProcess(const int &nodeid,string variablename,
                                         const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem){
    if(nodeid<1||nodeid>mesh.GetBulkMeshNodesNum()){
        MessagePrinter::PrintErrorTxt("nodeid="+to_string(nodeid)+" is invalid for NodeValuePostProcess");
        MessagePrinter::AsFem_Exit();
    }
    double nodevalue=0.0;
    int DofIndex,j;
    DofIndex=dofHandler.GetDofIDviaDofName(variablename);
    if(DofIndex<1){
        MessagePrinter::PrintErrorTxt("error detected in NodeValuePostProcess,"
                                      "we can not find DoF(name="+variablename+")"
                                      +", please check either your input file");
        MessagePrinter::AsFem_Exit();
    }

    j=dofHandler.GetBulkMeshIthNodeJthDofIndex(nodeid,DofIndex)-1;

    VecScatterCreateToAll(solutionSystem._Unew,&_scatteru,&_Useq);
    VecScatterBegin(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    VecGetValues(_Useq,1,&j,&nodevalue);
    // delete scatter
    VecScatterDestroy(&_scatteru);
    VecDestroy(&_Useq);
    return nodevalue;
}
