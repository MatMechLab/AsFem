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
//+++ Date   : 2020.07.01
//+++ Purpose: Implement general dofhandler for our bulk mesh
//+++          This class should be capable to manage DoFs, DoF maps...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"


BulkDofHandler::BulkDofHandler(){
    _nElmts=0;_nNodes=0;_nBulkElmts=0;
    _nDofsPerNode=0;
    _nDofs=0;_nActiveDofs=0;
    _HasDofMap=false;
    _HasSetDofName=false;
    _DofIDList.clear();
    _DofNameList.clear();
    _DofID2NameList.clear();
    _DofName2IDList.clear();

    vector<vector<int>> _NodeDofsMap;
    vector<vector<int>> _ElmtDofsMap;
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