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

#include "MateSystem/MaterialsContainer.h"

MaterialsContainer::MaterialsContainer(){
    m_BooleanMaterials.clear();
    m_ScalarMaterials.clear();
    m_VectorMaterials.clear();
    m_Rank2Materials.clear();
    m_Rank4Materials.clear();
}
void MaterialsContainer::clean(){
    m_BooleanMaterials.clear();
    m_ScalarMaterials.clear();
    m_VectorMaterials.clear();
    m_Rank2Materials.clear();
    m_Rank4Materials.clear();
}
MaterialsContainer::MaterialsContainer(const MaterialsContainer &a){
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
}
//*******************************************************
//*** for the access to each material
//*******************************************************
bool MaterialsContainer::BooleanMaterial(const string &matename)const{
    for(const auto &it:m_BooleanMaterials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("boolean material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return false;
}
double MaterialsContainer::ScalarMaterial(const string &matename)const{
    for(const auto &it:m_ScalarMaterials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("scalar material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return 0;
}
Vector3d MaterialsContainer::VectorMaterial(const string &matename)const{
    for(const auto &it:m_VectorMaterials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("vector material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return Vector3d(0);
}
Rank2Tensor MaterialsContainer::Rank2Material(const string &matename)const{
    for(const auto &it:m_Rank2Materials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("rank-2 tensor material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return Rank2Tensor(0.0);
}
Rank4Tensor MaterialsContainer::Rank4Material(const string &matename)const{
    for(const auto &it:m_Rank4Materials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("rank-4 tensor material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return Rank4Tensor(0.0);
}