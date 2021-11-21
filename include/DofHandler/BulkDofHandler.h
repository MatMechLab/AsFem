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
//+++ Date   : 2020.07.01
//+++ Update : 2021.11.13 @ Yang Bai
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

/**
 * this class implement dof management for the bulk element of AsFem.
 */
class BulkDofHandler{
public:
    /**
     * Constructor
     */
    BulkDofHandler();

    /**
     * add the dofs name from string vector
     * @param namelist the string vector for the name of dofs
     */
    void AddDofNameFromStrVec(vector<string> &namelist);

    /**
     * create the dofs map for bulk elements
     * @param mesh the mesh class
     * @param bcSystem the boundary condition system
     * @param elmtSystem the element system
     */
    void CreateBulkMeshDofsMap(const Mesh &mesh,BCSystem &bcSystem,ElmtSystem &elmtSystem);

    /**
     * get the i-th dof's name(here the dofs means the one defined in [dofs] block)
     * @param i the index of dofs in [dofs] block
     */
    inline string GetIthDofName(const int &i)const{return _DofNameList[i-1];}
    
    /**
     * get the i-t dof's ID(the id of dofs in [dofs] block)
     * @param i the dof index
     */
    inline int GetIthDofID(const int &i)const{return _DofIDList[i-1];}
   
    /**
     * get the nodes number of single bulk element
     */
    inline int GetNodesNumPerBulkElmt()const{return _nNodesPerBulkElmt;} 
    
    /**
     * get the maximum dofs number of each node
     */
    inline int GetMaxDofsNumPerNode()const{return _nDofsPerNode;}

    /**
     * get the maximum dofs number of each node
     */
    inline int GetDofsNumPerNode()const{return _nDofsPerNode;}
    
    /**
     * get the total dofs number of the system
     */
    inline int GetTotalDofsNum()const{return _nDofs;}
    
    /**
     * get the active dofs number of the system
     */
    inline int GetActiveDofsNum()const{return _nActiveDofs;}
    
    /**
     * get the maximum dofs number of single bulk element
     */
    inline int GetMaxDofsNumPerBulkElmt()const{return _nMaxDofsPerElmt;}
    
    /**
     * get te bulk elements number
     */
    inline int GetBulkMeshBulkElmtsNum()const{return _nBulkElmts;}
    
    /**
     * get the total elements number
     */
    inline int GetBulkMeshElmtsNum()const{return _nElmts;}

    /**
     * et the max non-zero entities of the row
     */
    inline int GetMaxRowNNZ()const{return _RowMaxNNZ;}
    
    /**
     * get the dof id by its name
     * @param dofname the string name of the inquery dof
     */
    int GetDofIDviaDofName(string dofname)const;

    /**
     * get the id index's vector from input dofs' name
     * @param namelist the input name string vector
     */
    vector<int> GetDofsIndexFromNameVec(vector<string> namelist)const;

    /**
     * get the i-th bulk element's dof index
     * @param e the bulk-element's ID
     * @param elDofs the int vector which stores the elemental dofs' ID
     * @param elDofsActiveFlag the active flags of the elemental dofs' ID, 1-> for active status,0-> for deactive status
     */
    inline void GetBulkMeshIthBulkElmtDofIndex(const int &e,vector<int> &elDofs,vector<double> &elDofsActiveFlag)const{
        for(int i=0;i<static_cast<int>(_BulkElmtDofsMap[e-1].size());i++){
            elDofs[i]=_BulkElmtDofsMap[e-1][i];
            elDofsActiveFlag[i]=_BulkElmtDofFlag[e-1][i];
        }
    }

    /**
     * get the i-th bulk element's dof index
     * @param e the bulk-element's ID
     * @param elDofs the int vector which stores the elemental dofs' ID
     */
    inline void GetBulkMeshIthBulkElmtDofIndex(const int &e,vector<int> &elDofs)const{
        for(int i=0;i<static_cast<int>(_BulkElmtDofsMap[e-1].size());i++){
            elDofs[i]=_BulkElmtDofsMap[e-1][i];
        }
    }

    /**
     * get the i-th bulk element's dof index and active flags, where the id index starts from 0!
     * @param e the bulk element's id
     * @param elDofs the elemental dofs' id list, start from 0
     * @param elDofsActiveFlag the elemental dofs' flag list
     */
    inline void GetBulkMeshIthBulkElmtDofIndex0(const int &e,vector<int> &elDofs,vector<double> &elDofsActiveFlag)const{
        for(int i=0;i<static_cast<int>(_BulkElmtDofsMap[e-1].size());i++){
            elDofs[i]=_BulkElmtDofsMap[e-1][i]-1;
            elDofsActiveFlag[i]=_BulkElmtDofFlag[e-1][i];
        }
    }

    /**
     * get the i-th bulk element's dof index and active flags, where the id index starts from 0!
     * @param e the bulk element's id
     * @param elDofs the elemental dofs' id list, start from 0
     */
    inline void GetBulkMeshIthBulkElmtDofIndex0(const int &e,vector<int> &elDofs)const{
        for(int i=0;i<static_cast<int>(_BulkElmtDofsMap[e-1].size());i++){
            elDofs[i]=_BulkElmtDofsMap[e-1][i]-1;
        }
    }

    /**
     * get the i-th bulk element's dofs id, which starts from 0!!!
     * @param e the bulk element's id
     * @param elDofs integer array's pointer
     */
    inline void GetBulkMeshIthBulkElmtDofIndex0(const int &e,int *elDofs)const{
        for(int i=0;i<static_cast<int>(_BulkElmtDofsMap[e-1].size());i++){
            elDofs[i]=_BulkElmtDofsMap[e-1][i]-1;
        }
    }
    
    /**
     * get the i-th bulk element's dofs id, which starts from 1!!!
     * @param e the bulk element's id
     * @param elDofs integer array's pointer
     */
    inline void GetBulkMeshIthBulkElmtDofIndex(const int &e,int *elDofs)const{
        for(int i=0;i<static_cast<int>(_BulkElmtDofsMap[e-1].size());i++){
            elDofs[i]=_BulkElmtDofsMap[e-1][i];
        }
    }

    /**
     * get i-th bulk element dofs number
     * @param e the element index, start from 1
     */
    inline int GetBulkMeshIthBulkElmtDofsNum(const int &e)const{
        return static_cast<int>(_BulkElmtDofsMap[e-1].size());
    }

    /**
     * get the i-th bulk element's j-th sub-element's dofs ID, start from 1 !!!
     * @param i bulk element id
     * @param j sub-element id
     */
    inline vector<int> GetBulkMeshIthBulkElmtJthSubElmtDofIndex(const int &i,const int &j)const{
        return _BulkElmtLocalDofIndex[i-1][j-1];
    }

    /**
     * get the i-th bulk element's j-th sub-element's material index
     * @param i bulk element id
     * @param j sub element id
     */
    inline int GetBulkMeshIthBulkElmtJthSubElmtMateIndex(const int &i,const int &j)const{
        return _BulkElmtElmtMateIndexList[i-1][j-1];
    }

    /**
     * get bulk mesh's i-th node's j-th coordinate
     * @param i node id, start from 1
     * @param j coordinate component, start from 1
     */
    inline int GetBulkMeshIthNodeJthDofIndex(const int &i,const int &j)const{
        return _NodeDofsMap[i-1][j-1];
    }

    /**
     * get bulk mesh's i-th node's j-th coordinate, start from 0
     * @param i node id, start from 0
     * @param j coordinate component, start from 0
     */
    inline int GetBulkMeshIthNodeJthDofIndex0(const int &i,const int &j)const{
        return _NodeDofsMap[i-1][j-1]-1;
    }

    //*********************************************
    //*** for some basic check functions
    //*********************************************
    /**
     * check wether the input string is a valid name for one of the [dofs] list
     * @param dofname input string for the name of inquery dof
     */
    bool IsValidDofName(string dofname)const;
    /**
     * check wether the namelist is a valid string vector for the [dofs] list
     * @param namelist the string vector for the inquery dofs
     */
    bool IsValidDofNameVec(vector<string> namelist)const;

    //*********************************************
    //*** for some basic getting functions
    //*********************************************
    /**
     * get the element type-material type pair of i-th bulk element
     * @param e the element id
     */
    inline vector<pair<ElmtType,MateType>> GetBulkMeshIthBulkElmtElmtMateTypePair(const int &e)const{
        return _BulkElmtElmtMateTypePairList[e-1];
    }

    /**
     * get the elmt type of i-th bulk element's j-th sub-element
     * @param i i-th bulk element
     * @param j j-th sub element
     */
    inline ElmtType GetBulkMeshIthBulkElmtJthSubElmtElmtType(const int &i,const int &j)const{
        return _BulkElmtElmtMateTypePairList[i-1][j-1].first;
    }

    /**
     * get the element type vector of i-th bulk element
     * @param e the element id
     */
    inline vector<ElmtType> GetBulkMeshIthBulkElmtElmtTypeVec(const int &e)const{
        vector<ElmtType> temp;temp.clear();
        for(auto it:_BulkElmtElmtMateTypePairList[e-1]){
            temp.push_back(it.first);
        }
        return temp;
    }

    /**
     * get the material type of i-th bulk element's j-th sub-element
     * @param i bulk element id
     * @param j sub element id
     */
    inline MateType GetBulkMeshIthBulkElmtJthSubElmtMateType(const int &i,const int &j)const{
        return _BulkElmtElmtMateTypePairList[i-1][j-1].second;
    }
    
    /**
     * get the material type vector of i-th bulk element
     * @param e bulk element id
     */
    inline vector<MateType> GetBulkMeshIthBulkElmtMateTypeVec(const int &e)const{
        vector<MateType> temp;temp.clear();
        for(auto it:_BulkElmtElmtMateTypePairList[e-1]){
            temp.push_back(it.second);
        }
        return temp;
    }

    /**
     * print out the information of the bulk dofs system
     */
    void PrintBulkDofInfo()const;
    /**
     * print out the details of the bulk dofs system 
     */
    void PrintBulkDofDetailInfo()const;

protected:
    //*************************************************
    //*** for basic dof information
    //*************************************************
    int _nElmts,_nBulkElmts;
    int _nDofsPerNode,_nMaxDofsPerNode;
    int _nDofs,_nActiveDofs;
    int _nNodesPerBulkElmt,_nNodes;
    int _nMaxDim,_nMinDim;
    int _nMaxDofsPerElmt;
    bool _HasDofMap,_HasSetDofName;
    vector<int>              _DofIDList;
    vector<string>           _DofNameList;
    vector<pair<int,string>> _DofID2NameList;
    vector<pair<string,int>> _DofName2IDList;

    vector<vector<int>> _NodeDofsMap;
    vector<vector<double>> _NodalDofFlag,_BulkElmtDofFlag;
    vector<vector<int>> _BulkElmtDofsMap;

    vector<vector<pair<ElmtType,MateType>>> _BulkElmtElmtMateTypePairList;
    vector<vector<int>> _BulkElmtElmtMateIndexList; 
    vector<vector<vector<int>>> _BulkElmtLocalDofIndex;

    // for the length of non-zero element per row
    vector<int> _RowNNZ;
    int _RowMaxNNZ; // the max non-zero elements of all the rows

};
