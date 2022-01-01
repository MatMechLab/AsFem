//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
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

#include "MateSystem/MateNameDefine.h"

#include "Utils/MessagePrinter.h"

using namespace std;

/**
 * this class responsible for the solution space operation, i.e., store the solution vectors, store/access the material properties.
 */
class SolutionSystem{
public:
    SolutionSystem();

    void AddSolutionNameFromVec(vector<string> vec){
        _DofNameList=vec;_HasDofNameList=true;
    }
    void AddProjectionNameFromVec(vector<string> vec){
        _ProjectionNameList=vec;_HasProjNameList=true;_IsProjection=true;
    }
    void AddScalarMateProjectionNameFromVec(vector<string> vec){
        _ScalarMateProjectionNameList=vec;_HasScalarMateProjName=true;_IsProjection=true;
    }
    void AddVectorMateProjectionNameFromVec(vector<string> vec){
        _VectorMateProjctionNameList=vec;_HasVectorMateProjName=true;_IsProjection=true;
    }
    void AddRank2MateProjectionNameFromVec(vector<string> vec){
        _Rank2MateProjectionNameList=vec;_HasRank2MateProjName=true;_IsProjection=true;
    }
    void AddRank4MateProjectionNameFromVec(vector<string> vec){
        _Rank4MateProjectionNameList=vec;_HasRank4MateProjName=true;_IsProjection=true;
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
    inline int GetProjIDViaName(string projname)const{
        for(int i=0;i<GetProjNumPerNode();i++){
            if(_ProjectionNameList[i]==projname){
                return i+1;
            }
        }
        return -1;
    }
    string GetIthProjName(const int &i)const{return _ProjectionNameList[i-1];}
    vector<string> GetProjNameVec()const{return _ProjectionNameList;}

    inline int GetScalarMateProjNumPerNode()const{return _nScalarProjPerNode;}
    inline int GetScalarMateIDViaName(string scalarname)const{
        for(int i=0;i<GetScalarMateProjNumPerNode();i++){
            if(_ScalarMateProjectionNameList[i]==scalarname){
                return i+1;
            }
        }
        return -1;
    }
    inline string GetIthScalarMateName(const int &i)const{
        return _ScalarMateProjectionNameList[i-1];
    }
    inline vector<string> GetScalarMateNameVec()const{return _ScalarMateProjectionNameList;}
    //****************************************
    inline int GetVectorMateProjNumPerNode()const{return _nVectorProjPerNode;}
    inline int GetVectorMateIDViaName(string scalarname)const{
        for(int i=0;i<GetVectorMateProjNumPerNode();i++){
            if(_VectorMateProjctionNameList[i]==scalarname){
                return i+1;
            }
        }
        return -1;
    }
    inline string GetIthVectorMateName(const int &i)const{
        return _VectorMateProjctionNameList[i-1];
    }
    inline vector<string> GetVectorMateNameVec()const{return _VectorMateProjctionNameList;}
    //*******************************************
    inline int GetRank2MateProjNumPerNode()const{return _nRank2ProjPerNode;}
    inline int GetRank2MateIDViaName(string scalarname)const{
        for(int i=0;i<GetRank2MateProjNumPerNode();i++){
            if(_Rank2MateProjectionNameList[i]==scalarname){
                return i+1;
            }
        }
        return -1;
    }
    inline string GetIthRank2MateName(const int &i)const{
        return _Rank2MateProjectionNameList[i-1];
    }
    inline vector<string> GetRank2MateNameVec()const{return _Rank2MateProjectionNameList;}
    //*********************************************
    inline int GetRank4MateProjNumPerNode()const{return _nRank4ProjPerNode;}
    inline int GetRank4MateIDViaName(string scalarname)const{
        for(int i=0;i<GetRank4MateProjNumPerNode();i++){
            if(_Rank4MateProjectionNameList[i]==scalarname){
                return i+1;
            }
        }
        return -1;
    }
    inline string GetIthRank4MateName(const int &i)const{
        return _Rank4MateProjectionNameList[i-1];
    }
    inline vector<string> GetRank4MateNameVec()const{return _Rank4MateProjectionNameList;}
    //**********************************************
    bool IsProjection()const{return _IsProjection;}

    void UpdateMaterials();

    void PrintProjectionInfo()const;

    void ReleaseMem();

public:
    Vec _Unew,_U,_V;
    Vec _Utemp;
    Vec _Uold,_Vold;

    Vec _Proj;// this is used for projection quantities in each element
    Vec _ProjScalarMate,_ProjVectorMate,_ProjRank2Mate,_ProjRank4Mate;// this is used for the material properties

    // store all the material properties on each gauss point
    // this is different from the ProjMaterials, the ProjMaterials
    // store the nodal material properties(projected from gauss point to nodal point for output)
    // while this one stores all the properties on each gauss point !!!
    vector<ScalarMateType> _ScalarMaterials,_ScalarMaterialsOld;
    vector<VectorMateType> _VectorMaterials,_VectorMaterialsOld;
    vector<Rank2MateType> _Rank2TensorMaterials,_Rank2TensorMaterialsOld;
    vector<Rank4MateType> _Rank4TensorMaterials,_Rank4TensorMaterialsOld;


private:
    bool _IsInit=false,_IsProjection=false;
    vector<string> _DofNameList;
    vector<string> _ProjectionNameList;
    vector<string> _ScalarMateProjectionNameList,_VectorMateProjctionNameList;
    vector<string> _Rank2MateProjectionNameList,_Rank4MateProjectionNameList;

    int _nHistPerGPoint,_nProjPerNode;
    int _nScalarProjPerNode,_nVectorProjPerNode,_nRank2ProjPerNode,_nRank4ProjPerNode;
    int _nGPointsPerBulkElmt;
    int _nDofs,_nNodes,_nElmts;
    bool _HasProjNameList=false;
    bool _HasScalarMateProjName=false;
    bool _HasVectorMateProjName=false;
    bool _HasRank2MateProjName=false;
    bool _HasRank4MateProjName=false;
    bool _HasDofNameList=false;

};
