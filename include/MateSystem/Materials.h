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
//+++ Date   : 2021.08.12
//+++ Purpose: Implement the materials class for AsFem, this class
//+++          contains the scalar, vector, rank-2, and rank4 type
//+++          materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include "Utils/MessagePrinter.h"
#include "MateSystem/MateNameDefine.h"

#include "Utils/Vector3d.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

using namespace std;

/**
 * This class offers the general access to different material, i.e., 
 * scalar materials, vector materials, rank-2 materials... via 
 * the name of the material
 */
class Materials{
public:
    /**
     * Constructor for materials, all the variables will be initialized!
     */
    Materials();
    /**
     * Constructor with a material class
     * @param newmate the right-hand side material class
     */
    Materials(const Materials &newmate);

    /**
     * The '=' operator between two materials,
     * left-hand side materials will be overwrite by the right-hand side one( both the size and material name!!!)
     * @param newmate the right-hand side material name
     */
    inline Materials& operator=(const Materials &newmate){
        Clean();
        for(const auto &it:newmate._ScalarMaterials){
            _ScalarMaterials[it.first]=it.second;
        }
        for(const auto &it:newmate._VectorMaterials){
            _VectorMaterials[it.first]=it.second;
        }
        for(const auto &it:newmate._Rank2Materials){
            _Rank2Materials[it.first]=it.second;
        }
        for(const auto &it:newmate._Rank4Materials){
            _Rank4Materials[it.first]=it.second;
        }
        return *this;
    }

    /**
     * Get the reference to scalar materials
     */
    ScalarMateType& GetScalarMatePtr(){return _ScalarMaterials;}

    /**
     * Get the copy of scalar materials 
     */
    ScalarMateType  GetScalarMate()const{return _ScalarMaterials;}
    
    /**
     * Get the reference to vector materials 
     */
    VectorMateType& GetVectorMatePtr(){return _VectorMaterials;}
    
    /**
     * Get the copy of vector material
     */
    VectorMateType GetVectorMate()const{return _VectorMaterials;}

    /**
     * Get the reference to rank-2 materials
     */
    Rank2MateType&  GetRank2MatePtr(){return _Rank2Materials;}
   
    /**
     * Get the copy of rank-2 materials
     */ 
    Rank2MateType GetRank2Mate()const{return _Rank2Materials;}

    /**
     * Get the reference to rank-4 materials 
     */
    Rank4MateType&  GetRank4MatePtr(){return _Rank4Materials;}

    /**
     * Get the copy of rank-4 materials
     */
    Rank4MateType  GetRank4Mate()const{return _Rank4Materials;}

    /**
     * This function will clean all the materials
     */
    void Clean();
    
    /**
     * Get the reference of the scalar type materials
     * @param matename the name of the material, if it is not there, AsFem will create one for you
     */
    double& ScalarMaterials(string matename);
    
    /**
     * Get the scalar material's value via its material name
     * @param matename the name of the material
     */
    double ScalarMaterials(string matename) const;
    
    /**
     * Get the refence of the vector type material
     * @param matename the material property's name
     */
    Vector3d& VectorMaterials(string matename);
    
    /**
     * Get the vector value of the vector type material, its always 3x1 size
     * @param matename the material property's name
     */
    Vector3d VectorMaterials(string matename) const;
    
    /**
     * Get the refence of the rank-2 type material
     * @param matename the material property's name
     */
    RankTwoTensor& Rank2Materials(string matename);
    
    /**
     * Get the rank-2 tensor value of the rank-2 tensor material, its always 3x3 size
     * @param matename the material property's name
     */
    RankTwoTensor  Rank2Materials(string matename) const;
    
    /**
     * Get the refence of the rank-4 type material
     * @param matename the material property's name
     */
    RankFourTensor& Rank4Materials(string matename);
    
    /**
     * Get the rank-4 value of the rank-4 type material, its always 3x3x3x3 size
     * @param matename the material property's name
     */
    RankFourTensor  Rank4Materials(string matename) const;

private:
    ScalarMateType _ScalarMaterials;/** the scalar type materials */
    VectorMateType _VectorMaterials;/** the vector type materials*/
    Rank2MateType  _Rank2Materials;/** the rank-2 type materials*/
    Rank4MateType  _Rank4Materials;/** the rank-4 type materials*/

    Vector3d       _VectorNull;/** the null(zero) vector*/
    RankTwoTensor  _Rank2Null; /** the null(zero) rank2 tensor*/
    RankFourTensor _Rank4Null; /** the null(zero) rank4 tensor*/

};
