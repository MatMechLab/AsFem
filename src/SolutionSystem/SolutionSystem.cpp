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
//+++ Date   : 2022.06.13
//+++ Purpose: defines the solution array for FEM analysis, stores
//+++          the necessary results.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SolutionSystem/SolutionSystem.h"

SolutionSystem::SolutionSystem(){
    m_Dofs=0;
    m_Allocated=false;
    m_BulkElmtsNum=0;
    m_QpointsNum=0;
}

void SolutionSystem::releaseMemory(){
    if(m_Allocated){
        m_Ucurrent.releaseMemory();
        m_Uold.releaseMemory();
        m_Uolder.releaseMemory();
        m_Utemp.releaseMemory();
        m_Ucopy.releaseMemory();

        m_V.releaseMemory();

        m_A.releaseMemory();

        // for local one
        for(auto &it:m_QpointsScalarMaterials_Local) it.clear();
        m_QpointsScalarMaterials_Local.clear();

        for(auto &it:m_QpointsVectorMaterials_Local) it.clear();
        m_QpointsVectorMaterials_Local.clear();

        for(auto &it:m_QpointsRank2Materials_Local) it.clear();
        m_QpointsRank2Materials_Local.clear();

        for(auto &it:m_QpointsRank4Materials_Local) it.clear();
        m_QpointsRank4Materials_Local.clear();
        // for total one
        for(auto &it:m_QpointsScalarMaterials_Total) it.clear();
        m_QpointsScalarMaterials_Total.clear();

        for(auto &it:m_QpointsVectorMaterials_Total) it.clear();
        m_QpointsVectorMaterials_Total.clear();

        for(auto &it:m_QpointsRank2Materials_Total) it.clear();
        m_QpointsRank2Materials_Total.clear();

        for(auto &it:m_QpointsRank4Materials_Total) it.clear();
        m_QpointsRank4Materials_Total.clear();

        //****************************************
        // for local one
        for(auto &it:m_QpointsScalarMaterialsOld_Local) it.clear();
        m_QpointsScalarMaterialsOld_Local.clear();

        for(auto &it:m_QpointsVectorMaterialsOld_Local) it.clear();
        m_QpointsVectorMaterialsOld_Local.clear();

        for(auto &it:m_QpointsRank2MaterialsOld_Local) it.clear();
        m_QpointsRank2MaterialsOld_Local.clear();

        for(auto &it:m_QpointsRank4MaterialsOld_Local) it.clear();
        m_QpointsRank4MaterialsOld_Local.clear();
        // for total one
        for(auto &it:m_QpointsScalarMaterialsOld_Total) it.clear();
        m_QpointsScalarMaterialsOld_Total.clear();

        for(auto &it:m_QpointsVectorMaterialsOld_Total) it.clear();
        m_QpointsVectorMaterialsOld_Total.clear();

        for(auto &it:m_QpointsRank2MaterialsOld_Total) it.clear();
        m_QpointsRank2MaterialsOld_Total.clear();

        for(auto &it:m_QpointsRank4MaterialsOld_Total) it.clear();
        m_QpointsRank4MaterialsOld_Total.clear();

        m_Allocated=false;
    }
}