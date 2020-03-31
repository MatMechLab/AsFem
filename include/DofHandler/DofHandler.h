//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_DOFHANDLER_H
#define ASFEM_DOFHANDLER_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

#include "petsc.h"

//***************************************
//*** For AsFem's own header file
//***************************************
#include "Mesh/Mesh.h"
#include "BCs/BCBlock.h"
#include "BCs/BCSystem.h"

using namespace std;


class Mesh;
class BCSystem;

class DofHandler{
public:
    DofHandler();
    DofHandler(const PetscInt &ndofspernode);
    //*********************************************************
    //*** add dofs' name to a list(to apply elmts, bc, ic or other stuff)
    //*********************************************************
    void AddNameToDofNameList(vector<string> dofnames);
    //*********************************************************
    //*** Create dof maps from mesh and boundary condition
    //*********************************************************
    void Init(Mesh &mesh);//do some basic allocation for related array
    void CreateDofMap(Mesh &mesh,BCSystem &bcSystem);// create dof maps for AsFem

    //*********************************************************
    //*** basic settings for dofhandler
    //*********************************************************
    void SetDofsNumPerNode(const PetscInt &ndofs){_nMaxDofsPerNode=ndofs;}
    void SetNodesNum(const PetscInt &nnodes){_nNodes=nnodes;}
    void SetElmtsNum(const PetscInt &nelmts){_nElmts=nelmts;}
    void SetBulkElmtsNum(const PetscInt &nelmts){_nBulkElmts=nelmts;}
    //*********************************************************
    //*** basic getting for dofhandler
    //*********************************************************
    inline PetscInt GetDofsNumPerNode()const{return _nMaxDofsPerNode;}
    inline PetscInt GetDofsNum()const{return _nDofs;}//=nNodes*nMaxDofsPerNode,but some domain may has less dofs per node than the max one!!!
    inline PetscInt GetDofsNumPerElmt()const{return _nMaxDofsPerElmt;}
    inline PetscInt GetMaxDofsNumPerElmt()const{return _nMaxDofsPerElmt;}
    inline PetscInt GetActiveDofsNum()const{return _nActiveDofs;}// the real dofs number we use
    inline PetscInt GetNodesNum()const{return _nNodes;}
    inline PetscInt GetElmtsNum()const{return _nElmts;}
    inline PetscInt GetBulkElmtsNum()const{return _nBulkElmts;}
    inline string GetIthDofNameFromList(const int &i)const{
        return _DofNameList[i-1];
    }
    //**** For NNZ in sparse matrix
    inline PetscInt GetMaxRowNNZ()const{
        return _RowMaxNNZ;
    }
    //************************
    //*** get dof index
    //************************
    vector<string> GetDofNameList()const{
        return _DofNameList;
    }
    inline PetscInt GetDofIndexViaName(string dofname)const{
        for(auto it:_DofNameToIDMap){
            if(it.first==dofname){
                return it.second;
            }
        }
        return 0;
    }
    //*** get dofs index from input string vector
    inline vector<PetscInt> GetDofsIndexFromNameVec(vector<string> namevec)const{
        vector<PetscInt> temp;
        temp.clear();
        for(auto it:namevec){
            temp.push_back(GetDofIndexViaName(it));
        }
        return temp;
    }
    //***
    bool IsValidDofName(string dofname)const{
        bool IsValid=false;
        for(auto it:_DofNameList){
            if(it==dofname){
                IsValid=true;
                break;
            }
        }
        return IsValid;
    }
    //***
    bool IsValidDofNameVec(vector<string> dofnamevec)const{
        bool IsValid=true;
        for(auto it:dofnamevec){
            IsValid=false;
            for(auto itj:_DofNameList){
                if(it==itj){
                    IsValid=true;
                    break;
                }
            }
            if(!IsValid){
                return false;
            }
        }
        return true;
    }
    //**************************************
    //*** for local dof index information
    //**************************************
    void GetIthBulkElmtDofIndex(const PetscInt &i,vector<PetscInt> &conn){
        for(int j=1;j<=(PetscInt)_ElmtalDofIndex[i-1].size();++j){
            conn[j-1]=_ElmtalDofIndex[i-1][j-1];
        }
    }
    
    void GetIthBulkElmtDofIndex0(const PetscInt &i,vector<PetscInt> &conn)const{
        for(int j=1;j<=(PetscInt)_ElmtalDofIndex[i-1].size();++j){
            conn[j-1]=_ElmtalDofIndex[i-1][j-1]-1;
        }
    }
    void GetIthBulkElmtDofIndex0(const PetscInt &i,vector<PetscInt> &conn,vector<PetscReal> dofactive)const{
        for(int j=1;j<=(PetscInt)_ElmtalDofIndex[i-1].size();++j){
            conn[j-1]=_ElmtalDofIndex[i-1][j-1]-1;
            dofactive[j-1]=_ElmtalDofActiveFlag[i-1][j-1];
        }
    }
    void GetIthBulkElmtDofIndex0(const PetscInt &i,PetscInt (&conn)[27])const{
        for(int j=1;j<=(PetscInt)_ElmtalDofIndex[i-1].size();++j){
            conn[j-1]=_ElmtalDofIndex[i-1][j-1]-1;
        }
    }
    inline PetscInt GetIthBulkElmtDofsNum(const PetscInt &i)const{
        return (PetscInt)_ElmtalDofIndex[i-1].size();
    }
    inline PetscInt GetIthNodeJthDofIndex(const PetscInt &i,const PetscInt &j)const{
        return _NodalDofFlag[i-1][j-1];
    }

public:
    //**********************************
    //*** print information
    //**********************************
    void PrintDofInfo() const;

private:
    PetscInt _nMaxDim,_nMinDim;
    PetscInt _nMaxDofsPerNode;
    PetscInt _nMaxDofsPerElmt,_nNodesPerBulkElmt;// for bulk elmt
    PetscInt _nDofsPerSurfaceElmt,_nDofsPerLineElmt;
    PetscInt _nDofs,_nActiveDofs;
    PetscInt _nBulkElmts,_nElmts,_nNodes;
    vector<string> _DofNameList;
    vector<string> _DofIndexList;
    vector<pair<string,PetscInt>> _DofNameToIDMap;
    vector<pair<PetscInt,string>> _DofIDToNameMap;
private:
    bool _IsDofMapCreated;
    //**********************************************
    //* for nodal and elemental dof flag
    //* for nodal doflag, each node has the max dofs pernode as the inital condition
    //* then, if node i's j-th node is dirichlet bc, then NodalDofFlag[i][j]=-1
    //* otherwise NodalDofFlag[i][j]!=-1, in the meantime, the value of NodalDofFlag[i][j] is the
    //* global dof index number(if it is dirichlet, then is -1)
    //* consequently, the elmtaldofflag just store the related dofs of related nodes in current elmts
    vector<vector<PetscInt>> _NodalDofFlag;
    vector<vector<PetscReal>> _NodalDofActiveFlag,_ElmtalDofActiveFlag;
    vector<vector<PetscInt>> _ElmtalDofIndex;
    vector<PetscInt> _RowNNZ;
    PetscInt _RowMaxNNZ;
};

#endif