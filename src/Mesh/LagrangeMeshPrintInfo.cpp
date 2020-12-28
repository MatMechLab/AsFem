//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.30
//+++ Purpose: Here we offer a printer to show the summarized info
//+++          for our mesh generation(either generated or imported)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/LagrangeMesh.h"

void LagrangeMesh::PrintBulkMeshInfo()const{
    char buff[70];
    // MessagePrinter::PrintDashLine();
    MessagePrinter::PrintNormalTxt("Mesh information summary:");
    
    snprintf(buff,70,"  nodes=%9d, nodesperbulkelmt=%3d, dim=%1d(max)-->%1d(min)",GetBulkMeshNodesNum(),GetBulkMeshNodesNumPerBulkElmt(),GetBulkMeshDim(),GetBulkMeshMinDim());
    MessagePrinter::PrintNormalTxt(string(buff));


    snprintf(buff,70,"  elmts=%9d, bulk elmts=%9d, sub dim elmts=%6d",GetBulkMeshElmtsNum(),GetBulkMeshBulkElmtsNum(),GetBulkMeshElmtsNum()-GetBulkMeshBulkElmtsNum());
    MessagePrinter::PrintNormalTxt(string(buff));


    snprintf(buff,70,"  physical groups=%4d, bulk mesh type=%8s, bulk mesh order=%2d",
                    GetBulkMeshPhysicalGroupNum(),GetBulkMeshBulkElmtTypeName().c_str(),GetBulkMeshOrder());
    MessagePrinter::PrintNormalTxt(string(buff));

    MessagePrinter::PrintNormalTxt("  physical id       dim                phsical Name          elmts");

    int phyid,n,dim;
    string phyname;
    for(int i=0;i<GetBulkMeshPhysicalGroupNum();i++){
        phyid=GetBulkMeshIthPhysicalID(i+1);
        phyname=GetBulkMeshIthPhysicalName(i+1);
        n=GetBulkMeshElmtsNumViaPhysicalName(phyname);
        dim=GetBulkMeshDimViaPhyName(phyname);
        snprintf(buff,70,"  %6d            %2d %28s      %8d",phyid,dim,phyname.c_str(),n);
        MessagePrinter::PrintNormalTxt(string(buff));
    }
    MessagePrinter::PrintDashLine();
}
//*********************************
void LagrangeMesh::PrintBulkMeshInfoDetails()const{
    char buff[70];
    MessagePrinter::PrintDashLine();
    MessagePrinter::PrintNormalTxt("Mesh information summary:");
    
    snprintf(buff,70,"  nodes=%9d, elmts=%9d, nodesperbulkelmt=%3d",GetBulkMeshNodesNum(),GetBulkMeshElmtsNum(),GetBulkMeshNodesNumPerBulkElmt());
    MessagePrinter::PrintNormalTxt(string(buff));

    snprintf(buff,70,"  max dim=%2d, min dim=%2d, phygroup=%4d, meshtype=%6s, order=%1d",
                     GetBulkMeshDim(),GetBulkMeshMinDim(),GetBulkMeshPhysicalGroupNum(),GetBulkMeshBulkElmtTypeName().c_str(),GetBulkMeshOrder());
    MessagePrinter::PrintNormalTxt(string(buff));

    MessagePrinter::PrintNormalTxt("  physical id                    phsical Name             elmts");

    int phyid,n;
    string phyname;
    for(int i=0;i<GetBulkMeshPhysicalGroupNum();i++){
        phyid=GetBulkMeshIthPhysicalID(i+1);
        phyname=GetBulkMeshIthPhysicalName(i+1);
        n=GetBulkMeshElmtsNumViaPhysicalName(phyname);
        snprintf(buff,70,"  %6d      %30s          %8d",phyid,phyname.c_str(),n);
        MessagePrinter::PrintNormalTxt(string(buff));
    }

    MessagePrinter::PrintNormalTxt("  Physical group information (ID and element ID) ");
    for(const auto &it:_PhysicalName2ElmtIDsList){
        for(int e=0;e<static_cast<int>(it.second.size());e++){
            snprintf(buff,70,"  phyname=%25s, element id=%9d",it.first.c_str(),it.second[e]);
            MessagePrinter::PrintNormalTxt(string(buff));
        }
    }

    char shortbuff[10],middlebuff[21];
    string str;
    MessagePrinter::PrintNormalTxt("  element connectivity information(element id: node index)");
    for(int e=1;e<=GetBulkMeshElmtsNum();++e){
        str.clear();
        snprintf(middlebuff,21,"  elmt id=%9d",e);
        str+=middlebuff;
        for(PetscInt i=1;i<=GetBulkMeshIthElmtNodesNum(e);++i){
            snprintf(shortbuff,10,"%8d",GetBulkMeshIthElmtJthNodeID(e,i));
            str+=shortbuff;
        }
        MessagePrinter::PrintNormalTxt(str);
    }

    MessagePrinter::PrintNormalTxt("node coornidates (node id, x, y, and z)");
    for(int i=1;i<=GetBulkMeshNodesNum();++i){
        snprintf(buff,70,"  %9d:%13.4e,%13.4e,%13.4e",i,GetBulkMeshIthNodeJthCoord(i,1),GetBulkMeshIthNodeJthCoord(i,2),GetBulkMeshIthNodeJthCoord(i,3));
        MessagePrinter::PrintNormalTxt(string(buff));
    }

    MessagePrinter::PrintDashLine();
}