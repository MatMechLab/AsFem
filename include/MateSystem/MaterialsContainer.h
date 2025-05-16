//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.07.25
//+++ Purpose: Implement the materials class for AsFem, this class
//+++          contains the scalar, vector, rank-2, and rank4 type
//+++          materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/MessagePrinter.h"
#include "MateSystem/MaterialsName.h"


/**
 * This class implemente the storage and management of scalar, vector, and tensor materials
 * of one single point or single element
 */
class MaterialsContainer{
public:
    /**
     * constructor
     */
    MaterialsContainer();
    /**
     * Constructor with a material container class
     * @param a the right-hand side material container class
     */
    MaterialsContainer(const MaterialsContainer &a);
    /**
     * clean up the current material container, all the material info will be removed
     */
    void clean();

    /**
     * The '=' operator between two materials,
     * left-hand side materials will be overwrite by the right-hand side one( both the size and material name!!!)
     * @param newmate the right-hand side material name
     */
    inline MaterialsContainer& operator=(const MaterialsContainer &a){
        clean();
        for(const auto &it:a.m_BooleanMaterials){
            m_BooleanMaterials[it.first]=it.second;
        }
        for(const auto &it:a.m_ScalarMaterials){
            m_ScalarMaterials[it.first]=it.second;
        }
        for(const auto &it:a.m_VectorMaterials){
            m_VectorMaterials[it.first]=it.second;
        }
        for(const auto &it:a.m_Rank2Materials){
            m_Rank2Materials[it.first]=it.second;
        }
        for(const auto &it:a.m_Rank4Materials){
            m_Rank4Materials[it.first]=it.second;
        }
        return *this;
    }
    //**************************************************************
    //*** general gettings
    //**************************************************************
    /**
     * get the boolean material value
     * @param matename the string name of inquery material
     */
    bool BooleanMaterial(const string &matename)const;
    /**
     * get the boolean material value
     * @param matename the string name of inquery material
     */
    inline bool& BooleanMaterial(const string &matename){
        return m_BooleanMaterials[matename];
    }

    /**
     * get the scalar material value
     * @param matename the string name of inquery material
     */
    double ScalarMaterial(const string &matename)const;
    /**
     * get the scalar material value
     * @param matename the string name of inquery material
     */
    inline double& ScalarMaterial(const string &matename){
        return m_ScalarMaterials[matename];
    }

    /**
     * get the vector material value
     * @param matename the string name of inquery material
     */
    Vector3d VectorMaterial(const string &matename)const;
    /**
     * get the vector material value
     * @param matename the string name of inquery material
     */
    inline Vector3d& VectorMaterial(const string &matename){
        return m_VectorMaterials[matename];
    }

    /**
     * get the rank-2 tensor material value
     * @param matename the string name of inquery material
     */
    Rank2Tensor Rank2Material(const string &matename)const;
    /**
     * get the rank-2 tensor material value
     * @param matename the string name of inquery material
     */
    Rank2Tensor& Rank2Material(const string &matename){
        return m_Rank2Materials[matename];
    }

    /**
     * get the rank-4 tensor material value
     * @param matename the string name of inquery material
     */
    Rank4Tensor Rank4Material(const string &matename)const;
    /**
     * get the rank-4 tensor material value
     * @param matename the string name of inquery material
     */
    Rank4Tensor& Rank4Material(const string &matename){
        return m_Rank4Materials[matename];
    }

    /**
     * get the reference of boolean materials
     */
    inline BooleanMateType& getBooleanMaterialsRef(){
        return m_BooleanMaterials;
    }
    /**
     * get the copy of boolean materials
     */
    inline BooleanMateType getBooleanMaterialsCopy()const{
        return m_BooleanMaterials;
    }

    /**
     * get the reference of scalar materials
     */
    inline ScalarMateType& getScalarMaterialsRef(){
        return m_ScalarMaterials;
    }
    /**
     * get the copy of scalar materials
     */
    inline ScalarMateType getScalarMaterialsCopy()const{
        return m_ScalarMaterials;
    }

    /**
     * get the reference of vector materials
     */
    inline VectorMateType& getVectorMaterialsRef(){
        return m_VectorMaterials;
    }
    /**
     * get the copy of vector materials
     */
    inline VectorMateType getVectorMaterialsCopy()const{
        return m_VectorMaterials;
    }

    /**
     * get the reference of rank-2 materials
     */
    inline Rank2MateType& getRank2MaterialsRef(){
        return m_Rank2Materials;
    }
    /**
     * get the copy of rank-2 materials
     */
    inline Rank2MateType getRank2MaterialsCopy()const{
        return m_Rank2Materials;
    }

    /**
     * get the reference of rank-4 materials
     */
    inline Rank4MateType& getRank4MaterialsRef(){
        return m_Rank4Materials;
    }
    /**
     * get the copy of rank-4 materials
     */
    inline Rank4MateType getRank4MaterialsCopy()const{
        return m_Rank4Materials;
    }

    /**
     * get the number of boolean materials
     */
    inline int getBooleanMaterialsNum()const{
        return static_cast<int>(m_BooleanMaterials.size());
    }
    /**
     * get the number of scalar materials
     */
    inline int getScalarMaterialsNum()const{
        return static_cast<int>(m_ScalarMaterials.size());
    }
    /**
     * get the number of vector materials
     */
    inline int getVectorMaterialsNum()const{
        return static_cast<int>(m_VectorMaterials.size());
    }
    /**
     * get the number of rank-2 tensor materials
     */
    inline int getRank2MaterialsNum()const{
        return static_cast<int>(m_Rank2Materials.size());
    }
    /**
     * get the number of rank-4 materials
     */
    inline int getRank4MaterialsNum()const{
        return static_cast<int>(m_Rank4Materials.size());
    }
    

private:
    BooleanMateType m_BooleanMaterials;/**< boolean materials */
    ScalarMateType m_ScalarMaterials;/**< scalar materials */
    VectorMateType m_VectorMaterials;/**< vector materials */
    Rank2MateType  m_Rank2Materials;/**< rank-2 materials */
    Rank4MateType  m_Rank4Materials;/**< rank-4 materials */

};
