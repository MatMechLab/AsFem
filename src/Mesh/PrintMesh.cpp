//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.08
//+++ Purpose: print out the bulk mesh info
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/BulkMesh.h"

void BulkMesh::printBulkMeshInfo()const{
    char buff[69];
    string str;
    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Mesh information summary");
    snprintf(buff,69,"  nodes=%8d, nodes per bulk elmt=%3d, max dim=%2d, min dim=%2d",getBulkMeshNodesNum(),
                                                                                    getBulkMeshNodesNumPerBulkElmt(),
                                                                                    getBulkMeshMaxDim(),
                                                                                    getBulkMeshMinDim());
    str=string(buff);
    MessagePrinter::printNormalTxt(str);

    snprintf(buff,69,"  elmts=%8d, bulk=%8d, surf=%8d, line=%8d",getTotalElmtsNum(),
                                                               getBulkMeshBulkElmtsNum(),
                                                               getBulkMeshSurfaceElmtsNum(),
                                                               getBulkMeshLineElmtsNum());
    str=buff;
    MessagePrinter::printNormalTxt(str);
    MessagePrinter::printNormalTxt("  bulk mesh type is "+getBulkMeshBulkElmtMeshTypeName()
                                  +", mesh order= "+to_string(getBulkMeshBulkElmtOrder()));

    snprintf(buff,69,"  total physical group=%5d, nodeset physical group=%5d",getBulkMeshPhyGroupNum(),getBulkMeshNodalPhyGroupNum());
    str=buff;
    MessagePrinter::printNormalTxt(str);

    MessagePrinter::printDashLine();

    MessagePrinter::printNormalTxt("  phy id              phy name       dim       nodes/elmt       elmts");
    for(int i=1;i<=getBulkMeshPhyGroupNum();i++){
        snprintf(buff,69,"%6d   %20s        %2d           %3d       %8d",getBulkMeshIthPhyGroupPhyID(i),
                                                 getBulkMeshIthPhyGroupPhyName(i).c_str(),
                                                 getBulkMeshIthPhyGroupDim(i),
                                                 getBulkMeshIthPhyGroupElmtNodesNum(i),
                                                 getBulkMeshIthPhyGroupElmtsNum(i));
        str=buff;
        MessagePrinter::printNormalTxt(str);
    }

    if(getBulkMeshNodalPhyGroupNum()){
        MessagePrinter::printDashLine();
        MessagePrinter::printNormalTxt("  nodal phy id             phy name             nodes number    ");
    }
    for(int i=1;i<=getBulkMeshNodalPhyGroupNum();i++){
        snprintf(buff,69,"   %6d     %20s                 %8d",getBulkMeshIthNodalPhyGroupPhyID(i),
                                                               getBulkMeshIthNodalPhyGroupPhyName(i).c_str(),
                                                               getBulkMeshIthNodalPhyGroupNodesNum(i));
        str=buff;
        MessagePrinter::printNormalTxt(str);
    }
    MessagePrinter::printStars();

}
void BulkMesh::printBulkMeshDetailedInfo()const{

}