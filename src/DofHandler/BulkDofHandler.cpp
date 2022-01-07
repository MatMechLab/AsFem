//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.01
//+++ Purpose: Implement general dofhandler for our bulk mesh
//+++          This class should be capable to manage DoFs, DoF maps...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"


BulkDofHandler::BulkDofHandler(){
    _nElmts=0.0;_nNodes=0.0;_nBulkElmts=0;
    _nDofsPerNode=0;_nNodesPerBulkElmt=0;
    _nDofs=0;_nActiveDofs=0;
    _nNodesPerBulkElmt=0;
    _nMaxDim=0;_nMinDim=0;
    _nMaxDofsPerElmt=0;
    _HasDofMap=false;_HasSetDofName=false;

    _DofIDList.clear();
    _DofNameList.clear();
    _DofID2NameList.clear();
    _DofName2IDList.clear();

    _DofIDList.clear();
    _DofNameList.clear();
    _DofID2NameList.clear();
    _DofName2IDList.clear();

    _NodeDofsMap.clear();
    _BulkElmtDofsMap.clear();

    _BulkElmtElmtMateTypePairList.clear();
}

void BulkDofHandler::AddDofNameFromStrVec(vector<string> &namelist){
    vector<string> namelistcopy=namelist;
    sort(namelistcopy.begin(),namelistcopy.end());
    namelistcopy.erase(unique(namelistcopy.begin(),namelistcopy.end()),namelistcopy.end());

    if(namelistcopy.size()!=namelistcopy.size()){
        MessagePrinter::PrintErrorTxt("the input dof namelist is not unique! you have duplicate dofs name");
        MessagePrinter::AsFem_Exit();
    }
    _DofIDList.clear();
    _DofNameList.clear();
    _DofID2NameList.clear();
    _DofName2IDList.clear();

    int i=0;
    for(auto it:namelist){
        i+=1;
        _DofIDList.push_back(i);
        _DofNameList.push_back(it);
        _DofID2NameList.push_back(make_pair(i,it));
        _DofName2IDList.push_back(make_pair(it,i));
    }
    _nDofsPerNode=i;
}
//**************************************************
void BulkDofHandler::PrintBulkDofInfo()const{
    char buff[70];
    MessagePrinter::PrintNormalTxt("Degrees of freedom (DoFs) information summary");
    snprintf(buff,70,"  each node has %2d dofs (max), total dofs=%6d",GetDofsNumPerNode(),GetActiveDofsNum());
    MessagePrinter::PrintNormalTxt(string(buff));

    MessagePrinter::PrintNormalTxt("  DoFs id                                 DoFs name");
    char longbuff[67];
    for(auto &it:_DofID2NameList){
        snprintf(longbuff,67,"  %2d                       %20s",it.first,it.second.c_str());
        MessagePrinter::PrintNormalTxt(string(longbuff));
    }
    MessagePrinter::PrintDashLine();
}

//**************************************************
void BulkDofHandler::PrintBulkDofDetailInfo()const{
    PrintBulkDofInfo();
    string str;
    char buff[14];
    MessagePrinter::PrintNormalTxt("Degrees of freedom (DoFs) information summary in details");
    for(int e=0;e<_nBulkElmts;e++){
        str.clear();
        sprintf(buff,"e=%8d:",e+1);
        str+=string(buff);
        for(auto it:_BulkElmtDofsMap[e]){
            sprintf(buff,"%8d ",it);
            str+=string(buff);
        }
        MessagePrinter::PrintLongTxt(str);
    }
    MessagePrinter::PrintDashLine();
}
