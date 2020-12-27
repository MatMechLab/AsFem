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
#include <map>
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

    void CreateBulkDofsMap(const Mesh &mesh,BCSystem &bcSystem,ElmtSystem &elmtSystem);

    inline string GetIthDofName(const int &i)const{return _DofNameList[i-1];}
    inline int    GetIthDofID(const int &i)const{return _DofIDList[i-1];}
    inline int    GetDofsNumPerNode()const{return _nDofsPerNode;}
    inline int    GetDofsNum()const{return _nDofs;}
    inline int    GetActiveDofsNum()const{return _nActiveDofs;}
    inline int    GetMaxDofsNumPerBulkElmt()const{return _nMaxDofsPerElmt;}
    inline int    GetBulkElmtNums()const{return _nBulkElmts;}
    inline int    GetMaxRowNNZ()const{return _RowMaxNNZ;}
    int         GetDofIDviaDofName(string dofname)const;
    vector<int> GetDofsIndexFromNameVec(vector<string> namelist)const;

    inline void GetIthBulkElmtDofIndex(const int &e,vector<int> &elDofs,vector<double> &elDofsActiveFlag)const{
        for(int i=0;i<static_cast<int>(_ElmtDofsMap[e-1].size());i++){
            elDofs[i]=_ElmtDofsMap[e-1][i];
            elDofsActiveFlag[i]=_ElmtDofFlag[e-1][i];
        }
    }

    inline void GetIthBulkElmtDofIndex(const int &e,int (&elDofs)[27])const{
        for(int i=0;i<static_cast<int>(_ElmtDofsMap[e-1].size());i++){
            elDofs[i]=_ElmtDofsMap[e-1][i];
        }
    }

    inline int GetIthBulkElmtDofsNum(const int &e)const{
        return static_cast<int>(_ElmtDofsMap[e-1].size());
    }

    inline vector<int> GetIthBulkElmtJthKernelDofIndex(const int &i,const int &j)const{
        return _ElmtLocalDofIndex[i-1][j-1];
    }
    inline int GetIthBulkElmtJthKernelMateIndex(const int &i,const int &j)const{
        return _ElmtElmtMateIndexList[i-1][j-1];
    }

    inline int GetIthNodeJthDofIndex(const int &i,const int &j)const{
        return _NodeDofsMap[i-1][j-1];
    }

    //*********************************************
    //*** for some basic check functions
    //*********************************************
    bool IsValidDofName(string dofname)const;
    bool IsValidDofNameVec(vector<string> namelist)const;

    //*********************************************
    //*** for some basic getting functions
    //*********************************************
    inline ElmtType GetIthElmtJthKernelElmtType(const int &i,const int &j)const{
        return _ElmtElmtMateTypePairList[i-1][j-1].first;
    }
    inline vector<pair<ElmtType,MateType>> GetIthElmtElmtMateTypePair(const int &e)const{
        return _ElmtElmtMateTypePairList[e-1];
    }
    inline vector<ElmtType> GetIthElmtElmtTypeVec(const int &e)const{
        vector<ElmtType> temp;temp.clear();
        for(auto it:_ElmtElmtMateTypePairList[e-1]){
            temp.push_back(it.first);
        }
        return temp;
    }
    inline MateType GetIthElmtJthKernelMateType(const int &i,const int &j)const{
        return _ElmtElmtMateTypePairList[i-1][j-1].second;
    }
    inline vector<MateType> GetIthElmtMateTypeVec(const int &e)const{
        vector<MateType> temp;temp.clear();
        for(auto it:_ElmtElmtMateTypePairList[e-1]){
            temp.push_back(it.second);
        }
        return temp;
    }

    void PrintBulkDofInfo()const;
    void PrintBulkDofDetailInfo()const;

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
    vector<vector<double>> _NodalDofFlag,_ElmtDofFlag;
    vector<vector<int>> _ElmtDofsMap;

    vector<vector<pair<ElmtType,MateType>>> _ElmtElmtMateTypePairList;
    vector<vector<int>> _ElmtElmtMateIndexList; 
    vector<vector<vector<int>>> _ElmtLocalDofIndex;

    // for the length of non-zero element per row
    vector<int> _RowNNZ;
    int _RowMaxNNZ; // the max non-zero elements of all the rows

};