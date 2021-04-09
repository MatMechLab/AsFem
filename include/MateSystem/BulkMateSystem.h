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

//*** For all the materials classes
#include "MateSystem/DoubleWellFreeEnergyMaterial.h"
#include "MateSystem/LinearElasticMaterial.h"
#include "MateSystem/MieheFractureMaterial.h"

class BulkMateSystem: public DoubleWellFreeEnergyMaterial,
        public LinearElasticMaterial,
        public MieheFractureMaterial{
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
        return _Materials.ScalarMaterials;
    }
    inline ScalarMateType& GetScalarMateOldPtr(){
        return _MaterialsOld.ScalarMaterials;
    }
    inline VectorMateType& GetVectorMatePtr(){
        return _Materials.VectorMaterials;
    }
    inline VectorMateType& GetVectorMateOldPtr(){
        return _MaterialsOld.VectorMaterials;
    }
    inline Rank2MateType& GetRank2MatePtr(){
        return _Materials.Rank2Materials;
    }
    inline Rank2MateType& GetRank2MateOldPtr(){
        return _MaterialsOld.Rank2Materials;
    }
    inline Rank4MateType& GetRank4MatePtr(){
        return _Materials.Rank4Materials;
    }
    inline Rank4MateType& GetRank4MateOldPtr(){
        return _MaterialsOld.Rank4Materials;
    }
    inline Materials& GetMaterialsPtr(){
        return _Materials;
    }
    inline Materials& GetMaterialsOldPtr(){
        return _MaterialsOld;
    }
    //***************************************************************************
    //*** for each materials getting functions
    //***************************************************************************
    inline vector<string> GetScalarMateNameList()const{
        vector<string> temp;
        temp.clear();
        for(const auto &it:_Materials.ScalarMaterials)temp.push_back(it.first);
        return temp;
    }
    inline bool IsNameInScalarMate(string matename)const{
        for(const auto &it:_Materials.ScalarMaterials){
            if(it.first==matename) return true;
        }
        return false;
    }
    inline int GetScalarMateNums()const{
        return static_cast<int>(_Materials.ScalarMaterials.size());
    }
    //*** For vector materials
    inline vector<string> GetVectorMateNameList()const{
        vector<string> temp;
        temp.clear();
        for(const auto &it:_Materials.VectorMaterials) temp.push_back(it.first);
        return temp;
    }
    inline bool IsNameInVectorMate(string matename)const{
        for(const auto &it:_Materials.VectorMaterials){
            if(it.first==matename) return true;
        }
        return false;
    }
    inline int GetVectorMateNums()const{
        return static_cast<int>(_Materials.VectorMaterials.size());
    }
    //*** For rank-2 materials
    inline vector<string> GetRank2MateNameList()const{
        vector<string> temp;
        temp.clear();
        for(const auto &it:_Materials.Rank2Materials)temp.push_back(it.first);
        return temp;
    }
    inline bool IsNameInRank2Mate(string matename)const{
        for(const auto &it:_Materials.Rank2Materials){
            if(it.first==matename) return true;
        }
        return false;
    }
    inline int GetRank2MateNums()const{
        return static_cast<int>(_Materials.Rank2Materials.size());
    }
    //*** For rank-4 materials
    inline vector<string> GetRank4MateNameList()const{
        vector<string> temp;
        temp.clear();
        for(const auto &it:_Materials.Rank4Materials)temp.push_back(it.first);
        return temp;
    }
    inline bool IsNameInRank4Mate(string matename)const{
        for(const auto &it:_Materials.Rank4Materials){
            if(it.first==matename) return true;
        }
        return false;
    }
    inline int GetRank4MateNums()const{
        return static_cast<int>(_Materials.Rank4Materials.size());
    }
    //***************************************************************************
    //*** For AsFem's built-in materials and User-Defined-Materials (UMAT)    ***
    //***************************************************************************
    void RunBulkMateLibs(const MateType &imate,const int &mateindex,const int &nDim,
                    const double &t,const double &dt,const Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpUOld,
                    const vector<double> &gpUdot,const vector<double> &gpUdotOld,
                    const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUOld,
                    const vector<Vector3d> &gpGradUdot,const vector<Vector3d> &gpGradUdotOld);

    void InitBulkMateLibs(const MateType &imate,const int &mateindex,const int &nDim,const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpUdot,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUdot);



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
    Materials _Materials,_MaterialsOld;


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

    //******************************************************************************
    //*** for constant diffusivity material
    //******************************************************************************
    void ConstDiffusionMaterial(const int &nDim,const double &t,const double &dt,
                              const vector<double> &InputParams,
                              const Vector3d &gpCoord,
                              const vector<double> &gpU,const vector<double> &gpV,
                              const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                              vector<double> &gpHist,const vector<double> &gpHistOld);

    //******************************************************************************
    //*** for cahnhilliard material
    //******************************************************************************
    void CahnHilliardMaterial(const int &nDim,const double &t,const double &dt,
                                const vector<double> &InputParams,
                                const Vector3d &gpCoord,
                                const vector<double> &gpU,const vector<double> &gpV,
                                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                                vector<double> &gpHist,const vector<double> &gpHistOld);

    //******************************************************************************
    //*** for User-Defined-Materials (UMAT)
    //******************************************************************************
    void User1Material(const int &nDim,const double &t,const double &dt,
                               const vector<double> &InputParams,
                               const Vector3d &gpCoord,
                               const vector<double> &gpU,const vector<double> &gpV,
                               const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                               vector<double> &gpHist,const vector<double> &gpHistOld);
};