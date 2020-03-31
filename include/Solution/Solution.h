//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_SOLUTION_H
#define ASFEM_SOLUTION_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"



using namespace std;

class Solution{
public:
    Solution();
    void AddSolutionNameFromVec(vector<string> vec){
        _DofNameList=vec;_HasDofNameList=true;
    }
    void AddProjectionNameFromVec(vector<string> vec){
        _ProjectionNameList=vec;_HasProjNameList=true;
    }
    void InitSolution(const PetscInt &ndofs,const PetscInt &nelmts,const PetscInt &nnodes,const PetscInt &ngp);

    //**************************************
    //*** Basic settings
    //**************************************
    void SetHistNumPerGPoint(const PetscInt &nhist){_nHistPerGPoint=nhist;}
    void SetProjNumPerNode(const PetscInt &nproj){_nProjPerNode=nproj;}

    //**************************************
    //*** Basic get funs
    //**************************************
    inline PetscInt GetHistNumPerGPoint()const{
        return _nHistPerGPoint;
    }
    inline PetscInt GetProjNumPerNode()const{
        return _nProjPerNode;
    }
    string GetIthProjName(const int &i)const{
        return _ProjectionNameList[i-1];
    }
    vector<string> GetProjNameVec()const{
        vector<string> temp=_ProjectionNameList;
        return temp;
    }

public:
    Vec _Unew,_Utemp;
    Vec _Uold,_Uolder;
    Vec _V,_Vold,_Volder;
    Vec _dU;

    Vec _Hist,_HistOld,_Proj;

private:
    bool _IsInit=false;
    vector<string> _DofNameList;
    vector<string> _ProjectionNameList;

    PetscInt _nHistPerGPoint,_nProjPerNode;
    PetscInt _nGPointsPerBulkElmt;
    PetscInt _nDofs,_nNodes,_nElmts;
    bool _HasProjNameList=false;
    bool _HasDofNameList=false;

};


#endif // ASFEM_SOLUTION_H