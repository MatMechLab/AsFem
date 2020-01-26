#ifndef ASFEM_DOFHANDLER_H
#define ASFEM_DOFHANDLER_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>


//**********************************
//*** AsFem's own header file    ***
//**********************************
#include "Mesh/Mesh.h"
#include "BCs/BCSystem.h"
#include "ElmtSystem/ElmtSystem.h"

using namespace std;


class BCSystem;
class DofHandler
{
public:
    DofHandler();
    void Init(Mesh &mesh);
    void CreateDofMap(Mesh &mesh,BCSystem &bcSystem);
    void ModifyDofActiveMapViaBC(Mesh &mesh,BCSystem &bcSystem);

public:
    // settings for dofhandler
    void SetDofNameListFromVec(vector<string> &vec);
    void SetDofNum(int ndofs) {_nDofs=ndofs;}
    void SetNodesNum(int nnodes) {_nNodes=nnodes;}
    void SetNodesNumPerBulkElmt(int nnodes) {_nNodesPerBulkElmt=nnodes;}
    void SetElmtsNum(int nelmts) {_nElmts=nelmts;}
    void SetBulkElmtsNum(int nelmts) {_nBulkElmts=nelmts;}
    
    inline void SetIthNodeJthDofFlag(const int &i,const int &j,const int &flag) {
        _NodalDofFlag[i-1][j-1]=flag;
    }
    inline void SetIthElmtJthDofFlag(const int &e,const int &j,const int &flag){
        _ElmtalDofIndex[e-1][j-1]=flag;
    }


    void PrintDofInfo() const;
    void PrintDofMap() const;
    

public:
    // get dofhandler information
    inline int GetDim() const {return _nDim;}
    inline int GetDofsNum() const {return _nDofs;}
    inline int GetActiveDofsNum() const {return _nActiveDofs;}
    inline int GetDofsNumPerNode() const {return _nDofsPerNode;}
    inline int GetMaxDofsNumPerElmt() const {return _nMaxDofsPerElmt;}
    inline int GetBulkElmtsNum() const {return _nBulkElmts;}
    inline int GetElmtsNum() const {return _nElmts;}
    inline int GetNodesNum() const {return _nNodes;}
    inline int GetNodesNumPerBulkElmt() const {return _nNodesPerBulkElmt;}
    inline int GetIthNodeJthDofFlag(const int &i,const int &j) const {return _NodalDofFlag[i-1][j-1];}
    inline int GetIthNodeJthDofIndex(const int &i,const int &j) const {return _NodalDofFlag[i-1][j-1];}

    inline int GetIthElmtDofsNum(const int &e) const{
        return int(_ElmtalDofIndex[e-1].size());
    }
    inline int GetIthElmtJthDofIndex(const int &e,const int &j) const{
        return _ElmtalDofIndex[e-1][j-1];
    }
    inline vector<int> GetIthElmtDofIndex(const int &e) const{
        // must be bulk element
        vector<int> elDofs;
        int nDofsPerElmt=int(_ElmtalDofIndex[e-1].size());
        elDofs.resize(nDofsPerElmt,0);
        for(int i=0;i<nDofsPerElmt;++i){
            elDofs[i]=_ElmtalDofIndex[e-1][i];
        }
        return elDofs;
    }
    inline void GetIthElmtDofIndex(const int &e,vector<int> &elDofs) const {
        for(int i=0;i<int(_ElmtalDofIndex[e-1].size());++i){
            elDofs[i]=_ElmtalDofIndex[e-1][i];
        }
    }
    inline void GetIthElmtDofIndex(const int &e,vector<int> &elDofs,vector<double> &elDofsActive) const {
        for(int i=0;i<int(_ElmtalDofIndex[e-1].size());++i){
            elDofs[i]=_ElmtalDofIndex[e-1][i];
            elDofsActive[i]=_ElmtalDofActiveFlags[e-1][i];
        }
    }

    inline string GetIthDofNameFromList(int i) {return _DofNameList[i-1];}
    inline int GetDofNameListLenght() const {return (int)_DofNameList.size();}
    inline int GetDofIndexViaName(string name) const 
    {
        for(int i=1;i<=int(_DofNameList.size());i++)
        {
            if(_DofNameList[i-1]==name)
            {
                return i;
            }
        }
        return 0;
    }

    inline vector<int> GetDofsIndexFromNameVec(vector<string> namevec) const{
        vector<int> temp;
        temp.clear();
        for(unsigned int j=0;j<namevec.size();++j){
            for(int i=1;i<=int(_DofNameList.size());++i){
                if(_DofNameList[i-1]==namevec[j]){
                    temp.push_back(i);
                    break;
                }
            }
        }
        return temp;
    }

    inline vector<string> GetDofNameList() const 
    {
        vector<string> temp;
        temp.clear();
        for(unsigned int i=0;i<_DofNameList.size();i++)
        {
            temp.push_back(_DofNameList[i]);
        }
        return temp;
    }


    inline void GetIthElmtalDofIndex(const int &e,vector<int> &ind) const{
        ind.clear();
        ind=_ElmtalDofIndex[e-1];
    }


    void PrintDofsInfo() const;

private:
    int _nDofs,_nNodes,_nElmts,_nBulkElmts,_nDofsPerNode,_nDim;
    int _nNodesPerBulkElmt,_nMaxDofsPerElmt;
    int _nActiveDofs;
    vector<string> _DofNameList;
    vector<int> _DofIndexList;// the same length as _DofNameList, related one by one!!!
    vector<int> _NodalDofIndex;
    vector<vector<int>> _ElmtalDofIndex;
    vector<vector<double>> _NodalDofActiveFlag,_ElmtalDofActiveFlags;
    vector<vector<int>> _NodalDofFlag;
};

#endif