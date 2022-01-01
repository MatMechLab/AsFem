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
#include "MateSystem/Materials.h"

#include "Utils/Vector3d.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

//*** For all the materials classes
#include "MateSystem/ConstPoissonMaterial.h"
#include "MateSystem/ConstDiffusionMaterial.h"
#include "MateSystem/DoubleWellFreeEnergyMaterial.h"
#include "MateSystem/IdealSolutionFreeEnergyMaterial.h"
#include "MateSystem/LinearElasticMaterial.h"
#include "MateSystem/IncrementSmallStrainMaterial.h"
#include "MateSystem/NeoHookeanMaterial.h"
#include "MateSystem/SaintVenantMaterial.h"
#include "MateSystem/Plastic1DMaterial.h"
#include "MateSystem/J2PlasticityMaterial.h"
#include "MateSystem/MieheFractureMaterial.h"
#include "MateSystem/StressDecompositionMaterial.h"
#include "MateSystem/NeoHookeanPFFractureMaterial.h"
#include "MateSystem/KobayashiMaterial.h"
#include "MateSystem/DiffNeoHookeanMaterial.h"
#include "MateSystem/DiffusionFractureMaterial.h"
#include "MateSystem/LinearElasticCHMaterial.h"
//*** For user defined materials(UMAT)
#include "MateSystem/User1Material.h"
#include "MateSystem/User2Material.h"
#include "MateSystem/User3Material.h"
#include "MateSystem/User4Material.h"
#include "MateSystem/User5Material.h"
#include "MateSystem/User6Material.h"
#include "MateSystem/User7Material.h"
#include "MateSystem/User8Material.h"
#include "MateSystem/User9Material.h"
#include "MateSystem/User10Material.h"


/**
 * This class offers the calculation for different materials(constitutive laws), which will 
 * be used by the ElmtSystem.
 */
class BulkMateSystem: public ConstPoissonMaterial,
                      public ConstDiffusionMaterial,
                      public DoubleWellFreeEnergyMaterial,
                      public IdealSolutionFreeEnergyMaterial,
                      public LinearElasticMaterial,
                      public IncrementSmallStrainMaterial,
                      public NeoHookeanMaterial,
                      public SaintVenantMaterial,
                      public Plastic1DMaterial,
                      public J2PlasticityMaterial,
                      public MieheFractureMaterial,
                      public StressDecompositionMaterial,
                      public NeoHookeanPFFractureMaterial,
                      public KobayashiMaterial,
                      public DiffNeoHookeanMaterial,
                      public DiffusionFractureMaterial,
                      public LinearElasticCHMaterial,
                      public User1Material,
                      public User2Material,
                      public User3Material,
                      public User4Material,
                      public User5Material,
                      public User6Material,
                      public User7Material,
                      public User8Material,
                      public User9Material,
                      public User10Material{
public:
    /*
     * The constructor of BulkMateSystem.
     */
    BulkMateSystem();
    BulkMateSystem(const BulkMateSystem &newbulkmatesystem);
   
    /**
     * for right-hand assignment/copy operator
     */
    inline BulkMateSystem& operator=(const BulkMateSystem &newbulkmatesystem){
        _nBulkMateBlocks=newbulkmatesystem._nBulkMateBlocks;
        _BulkMateBlockList=newbulkmatesystem._BulkMateBlockList;
        _Materials=newbulkmatesystem._Materials;
        _MaterialsOld=newbulkmatesystem._MaterialsOld;
        return *this;
    }

    /**
     * Initialize the material system.
     */
    void InitBulkMateSystem();
    
    /**
     * add the [mate] subblock of the input file to the material block list.
     */
    void AddBulkMateBlock2List(MateBlock &mateblock);
    //********************************************
    //*** for some basic getting functions
    //********************************************
    /**
     * return the total [mate] sub block number.
     */
    inline int GetMateBlockNums()const{return _nBulkMateBlocks;}
    
    /**
     * get the copy of i-th [mate] sub block.
     * @param i the index of [mate] sub block, start from 1, not 0!!!
     */
    inline MateBlock GetIthMateBlock(const int &i)const{return _BulkMateBlockList[i-1];}
    
    /**
     * get the copy of [mate] block vector
     */
    inline vector<MateBlock> GetMateBlockVec()const{return _BulkMateBlockList;}


    /**
     * get the reference of scalar materials
     */
    inline ScalarMateType& GetScalarMatePtr(){
        return _Materials.GetScalarMatePtr();
    }
    /**
     * get the reference of old scalar materials
     */
    inline ScalarMateType& GetScalarMateOldPtr(){
        return _MaterialsOld.GetScalarMatePtr();
    }

    /**
     * get the reference of vector materials
     */
    inline VectorMateType& GetVectorMatePtr(){
        return _Materials.GetVectorMatePtr();
    }
    /**
     * get the reference of old vector materials
     */
    inline VectorMateType& GetVectorMateOldPtr(){
        return _MaterialsOld.GetVectorMatePtr();
    }
    
    /**
     * get the reference of rank-2 materials
     */
    inline Rank2MateType& GetRank2MatePtr(){
        return _Materials.GetRank2MatePtr();
    }
    /**
     * get the reference of old rank-2 materials
     */
    inline Rank2MateType& GetRank2MateOldPtr(){
        return _MaterialsOld.GetRank2MatePtr();
    }

    /**
     * get the reference of rank-4 materials
     */
    inline Rank4MateType& GetRank4MatePtr(){
        return _Materials.GetRank4MatePtr();
    }
    /**
     * get the reference of old rank-4 materials
     */
    inline Rank4MateType& GetRank4MateOldPtr(){
        return _MaterialsOld.GetRank4MatePtr();
    }

    /**
     * get the reference of materials class
     */
    inline Materials& GetMaterialsPtr(){
        return _Materials;
    }
    /**
     * get the reference of old materials class
     */
    inline Materials& GetMaterialsOldPtr(){
        return _MaterialsOld;
    }
    //***************************************************************************
    //*** for each materials getting functions
    //***************************************************************************
    /**
     * get the name list of scalar materials
     */
    inline vector<string> GetScalarMateNameList()const{
        vector<string> temp;
        temp.clear();
        for(const auto &it:_Materials.GetScalarMate())temp.push_back(it.first);
        return temp;
    }

    /**
     * check wether the name is a valid scalar material name
     * @param matename the name of the material to be checked
     */
    inline bool IsNameInScalarMate(string matename)const{
        for(const auto &it:_Materials.GetScalarMate()){
            if(it.first==matename) return true;
        }
        return false;
    }

    /**
     * get the total number of scalar materials
     */
    inline int GetScalarMateNums()const{
        return static_cast<int>(_Materials.GetScalarMate().size());
    }
    //*** For vector materials
    /**
     * get the name list of vector materials
     */
    inline vector<string> GetVectorMateNameList()const{
        vector<string> temp;
        temp.clear();
        for(const auto &it:_Materials.GetVectorMate()) temp.push_back(it.first);
        return temp;
    }
    /**
     * chekc wether the material name is a valid vector material's name
     * @param matename the name to be checked
     */
    inline bool IsNameInVectorMate(string matename)const{
        for(const auto &it:_Materials.GetVectorMate()){
            if(it.first==matename) return true;
        }
        return false;
    }
    /**
     * get the total number of vector materials
     */
    inline int GetVectorMateNums()const{
        return static_cast<int>(_Materials.GetVectorMate().size());
    }
    //*** For rank-2 materials
    /**
     * get the name list of rank-2 materials
     */
    inline vector<string> GetRank2MateNameList()const{
        vector<string> temp;
        temp.clear();
        for(const auto &it:_Materials.GetRank2Mate())temp.push_back(it.first);
        return temp;
    }
    /**
     * check wether the name is valid for a rank-2 material
     */
    inline bool IsNameInRank2Mate(string matename)const{
        for(const auto &it:_Materials.GetRank2Mate()){
            if(it.first==matename) return true;
        }
        return false;
    }
    /**
     * get the total number of rank-2 materials
     */
    inline int GetRank2MateNums()const{
        return static_cast<int>(_Materials.GetRank2Mate().size());
    }
    //*** For rank-4 materials
    /**
     * get the name list of rank-4 materials 
     */
    inline vector<string> GetRank4MateNameList()const{
        vector<string> temp;
        temp.clear();
        for(const auto &it:_Materials.GetRank4Mate())temp.push_back(it.first);
        return temp;
    }
    /**
     * check wether the name is a valid one for rank-4 materials
     */
    inline bool IsNameInRank4Mate(string matename)const{
        for(const auto &it:_Materials.GetRank4Mate()){
            if(it.first==matename) return true;
        }
        return false;
    }
    /**
     * get the total number of rank-4 materials
     */
    inline int GetRank4MateNums()const{
        return static_cast<int>(_Materials.GetRank4Mate().size());
    }


    //***************************************************************************
    //*** For AsFem's built-in materials and User-Defined-Materials (UMAT)    ***
    //***************************************************************************
    /**
     * calculate the material properties by running different bulk material calc, the bulk material calc should be registered within this function.
     * @param imate the type of the bulk material
     * @param mateindex the index of the material block, its the order you defined in your input file
     * @param elmtinfo the information of current local element
     * @param elmtsoln the solution and its gradient of current local element
     */
    void RunBulkMateLibs(const MateType &imate,const int &mateindex,
                         const LocalElmtInfo &elmtinfo,
                         const LocalElmtSolution &elmtsoln);
    
    
    /**
     * calculate the initial material properties by running different bulk material , the bulk material calc should be registered within this function.
     * @param imate the type of the bulk material
     * @param mateindex the index of the material block, its the order you defined in your input file
     * @param elmtinfo the information of current local element
     * @param elmtsoln the solution and its gradient of current local element
     */
    void InitBulkMateLibs(const MateType &imate,const int &mateindex,
                          const LocalElmtInfo &elmtinfo,
                          const LocalElmtSolution &elmtsoln);



    /**
     * print the information of bulk material system
     */
    void PrintBulkMateSystemInfo()const;


protected:
    int _nBulkMateBlocks;/**< the number of bulk [mate] block*/
    vector<MateBlock> _BulkMateBlockList;/**< the vector of bulk [mate] list*/

protected:
    Materials _Materials;/**< the materials class for current time step, it contains scalar,vector,rank-2,rank-4 materials */
    Materials _MaterialsOld;/**< the materials class from previous step */


};
