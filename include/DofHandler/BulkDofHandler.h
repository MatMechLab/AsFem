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

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>

#include "Utils/MessagePrinter.h"
#include "Mesh/Mesh.h"
#include "BCSystem/BCSystem.h"
#include "ElmtSystem/ElmtSystem.h"

using namespace std;

class Mesh;
class ElmtSystem;
class BCSystem;

class BulkDofHandler{
public:
    BulkDofHandler();

    void AddDofNameFromStrVec(vector<string> &namelist);

    void CreateBulkDofsMap(Mesh &mesh,BCSystem &bcSystem,ElmtSystem &elmtSystem);

    inline string GetIthDofName(const int &i)const{return _DofNameList[i-1];}
    inline int    GetIthDofID(const int &i)const{return _DofIDList[i-1];}
    inline int    GetDofsNumPerNode()const{return _nDofsPerNode;}
    inline int    GetDofsNum()const{return _nDofs;}
    inline int    GetActiveDofsNum()const{return _nActiveDofs;}
    
    int         GetDofIDviaDofName(string dofname)const;
    vector<int> GetDofsIndexFromNameVec(vector<string> namelist)const;

    //*********************************************
    //*** for some basic check functions
    //*********************************************
    bool IsValidDofName(string dofname)const;
    bool IsValidDofNameVec(vector<string> namelist)const;

    //*********************************************
    //*** for some basic getting functions
    //*********************************************
    inline pair<ElmtType,MateType> GetIthElmtElmtMateTypePair(const int &e)const{
        return _ElmtElmtMateTypePairList[e-1];
    }
    inline ElmtType GetIthElmtElmtType(const int &e)const{
        return _ElmtElmtMateTypePairList[e-1].first;
    }
    inline MateType GetIthElmtMateType(const int &e)const{
        return _ElmtElmtMateTypePairList[e-1].second;
    }

    void PrintDofInfo()const;

protected:
    //*************************************************
    //*** for basic dof information
    //*************************************************
    int _nElmts,_nNodes,_nBulkElmts;
    int _nDofsPerNode;
    int _nDofs,_nActiveDofs;
    int _nNodesPerBulkElmt;
    int _nMaxDim,_nMinDim;
    int _nMaxDofsPerElmt;
    bool _HasDofMap,_HasSetDofName;
    vector<int>              _DofIDList;
    vector<string>           _DofNameList;
    vector<pair<int,string>> _DofID2NameList;
    vector<pair<string,int>> _DofName2IDList;

    vector<vector<int>> _NodeDofsMap;
    vector<vector<double>> _NodalDofFlag;
    vector<vector<int>> _ElmtDofsMap;

    vector<pair<ElmtType,MateType>> _ElmtElmtMateTypePairList;

};