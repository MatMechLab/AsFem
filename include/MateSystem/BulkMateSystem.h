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
//+++ Date   : 2020.11.30
//+++ Purpose: Implement the materials system for our bulk element
//+++          in AsFem. It is different from another one, namely
//+++          the interface material system
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>


//**********************************
//*** For AsFem's own header file
//**********************************
#include "Utils/MessagePrinter.h"
#include "MateSystem/MateBlock.h"
#include "MateSystem/MateType.h"
#include "MateSystem/MateTypeDefine.h"

#include "Utils/Vector3d.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

class BulkMateSystem{
public:
    BulkMateSystem();
    void InitBulkMateSystem();
    void AddBulkMateBlock2List(MateBlock &mateblock);
    //********************************************
    //*** for some basic getting functions
    //********************************************
    inline int GetMateBlockNums()const{return _nBulkMateBlocks;}
    inline MateBlock GetIthMateBlock(const int &i)const{return _BulkMateBlockList[i-1];}
    inline vector<MateBlock> GetMateBlockVec()const{return _BulkMateBlockList;}


    inline ScalarMateType& GetScalarMatePtr(){
        return _ScalarMaterials;
    }
    inline VectorMateType& GetVectorMatePtr(){
        return _VectorMaterials;
    }
    inline Rank2MateType& GetRank2MatePtr(){
        return _Rank2Materials;
    }
    inline Rank4MateType& GetRank4MatePtr(){
        return _Rank4Materials;
    }
    //***************************************************************************
    //*** For AsFem's built-in materials and User-Defined-Materials (UMAT)    ***
    //***************************************************************************
    void RunBulkMateLibs(const MateType &imate,const int &mateindex,const int &nDim,
                    const double &t,const double &dt,
                    const Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                    vector<double> &gpHist,const vector<double> &gpHistOld);


    void PrintBulkMateSystemInfo()const;


protected:
    int _nBulkMateBlocks;
    vector<MateBlock> _BulkMateBlockList;

protected:
    //***************************************************************************
    //*** Here we define four types of materials which will be used by our    ***
    //*** element system, we use the 'map' to simplify the access of material ***
    //*** properties by its name
    //*** Material types:
    //***  1) Scalar type materials
    //***  2) Vector type materials, by default, we assume the vector3d type
    //***  3) Rank-2 tensor type materials, i.e. stresses and strains
    //***  4) Rank-4 tensor type materials, i.e. elasticity tensor
    //****************************************************************************
    ScalarMateType _ScalarMaterials;
    VectorMateType _VectorMaterials;
    Rank2MateType  _Rank2Materials;
    Rank4MateType  _Rank4Materials;


protected:
    //******************************************************************************
    //*** Here we list all of the built-in materials as well as UMAT
    //******************************************************************************
    void ConstPoissonMaterial(const int &nDim,const double &t,const double &dt,
                        const vector<double> &InputParams,
                        const Vector3d &gpCoord,
                        const vector<double> &gpU,const vector<double> &gpV,
                        const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                        vector<double> &gpHist,const vector<double> &gpHistOld);
};