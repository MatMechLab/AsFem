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
//+++ Date   : 2020.07.12
//+++ Purpose: define the solution system for AsFem, where all the 
//+++          results/material variables/projection quantities 
//+++          should be stored here by this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"

#include "Utils/MessagePrinter.h"

using namespace std;

class SolutionSystem{
public:
    SolutionSystem();

    void AddSolutionNameFromVec(vector<string> vec){
        _DofNameList=vec;_HasDofNameList=true;
    }
    void AddProjectionNameFromVec(vector<string> vec){
        _ProjectionNameList=vec;_HasProjNameList=true;
    }
    void InitSolution(const int &ndofs,const int &nelmts,const int &nnodes,const int &ngp);

    //**************************************
    //*** Basic settings
    //**************************************
    void SetHistNumPerGPoint(const int &nhist){_nHistPerGPoint=nhist;}
    void SetProjNumPerNode(const int &nproj){_nProjPerNode=nproj;}
    void SetProjectionStatus(bool flag){_IsProjection=flag;}

    //**************************************
    //*** Basic get funs
    //**************************************
    inline int GetHistNumPerGPoint()const{return _nHistPerGPoint;}
    inline int GetProjNumPerNode()const{return _nProjPerNode;}
    string GetIthProjName(const int &i)const{return _ProjectionNameList[i-1];}
    vector<string> GetProjNameVec()const{return _ProjectionNameList;}
    bool IsProjection()const{return _IsProjection;}


    void PrintProjectionInfo()const;

public:
    Vec _Unew,_Utemp;
    Vec _Uold,_Uolder;
    Vec _V,_Vold,_Volder;
    Vec _dU;

    Vec _Hist,_HistOld,_Proj;

private:
    bool _IsInit=false,_IsProjection=false;
    vector<string> _DofNameList;
    vector<string> _ProjectionNameList;

    int _nHistPerGPoint,_nProjPerNode;
    int _nGPointsPerBulkElmt;
    int _nDofs,_nNodes,_nElmts;
    bool _HasProjNameList=false;
    bool _HasDofNameList=false;

};