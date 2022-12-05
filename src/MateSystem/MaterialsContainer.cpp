//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
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
    m_boolean_materials.clear();
    m_scalar_materials.clear();
    m_vector_materials.clear();
    m_rank2_materials.clear();
    m_rank4_materials.clear();
}
void MaterialsContainer::clean(){
    m_boolean_materials.clear();
    m_scalar_materials.clear();
    m_vector_materials.clear();
    m_rank2_materials.clear();
    m_rank4_materials.clear();
}
MaterialsContainer::MaterialsContainer(const MaterialsContainer &a){
    clean();
    for(const auto &it:a.m_boolean_materials){
        m_boolean_materials[it.first]=it.second;
    }
    for(const auto &it:a.m_scalar_materials){
        m_scalar_materials[it.first]=it.second;
    }
    for(const auto &it:a.m_vector_materials){
        m_vector_materials[it.first]=it.second;
    }
    for(const auto &it:a.m_rank2_materials){
        m_rank2_materials[it.first]=it.second;
    }
    for(const auto &it:a.m_rank4_materials){
        m_rank4_materials[it.first]=it.second;
    }
}
//*******************************************************
//*** for the access to each material
//*******************************************************
bool MaterialsContainer::BooleanMaterial(const string &matename)const{
    for(const auto &it:m_boolean_materials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("boolean material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return false;
}
double MaterialsContainer::ScalarMaterial(const string &matename)const{
    for(const auto &it:m_scalar_materials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("scalar material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return 0;
}
Vector3d MaterialsContainer::VectorMaterial(const string &matename)const{
    for(const auto &it:m_vector_materials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("vector material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return Vector3d(0);
}
Rank2Tensor MaterialsContainer::Rank2Material(const string &matename)const{
    for(const auto &it:m_rank2_materials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("rank-2 tensor material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return Rank2Tensor(0.0);
}
Rank4Tensor MaterialsContainer::Rank4Material(const string &matename)const{
    for(const auto &it:m_rank4_materials){
        if(it.first==matename){
            return it.second;
        }
    }
    MessagePrinter::printErrorTxt("rank-4 tensor material ("+matename+") is not defined in your materials");
    MessagePrinter::exitAsFem();
    return Rank4Tensor(0.0);
}