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
//+++ Date   : 2022.07.29
//+++ Purpose: update the solution array and materials array
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SolutionSystem/SolutionSystem.h"

void SolutionSystem::updateSolution(){
    m_Uolder=m_Uold;
    m_Uold=m_Ucurrent;
}

void SolutionSystem::updateMaterialsSolution(){
    /**
     * Here we only need to update the local material vectors
     */
    for(int e=1;e<=m_BulkElmtsNum_Local;e++){
        for(int j=1;j<=m_QpointsNum;j++){
            m_QpointsScalarMaterialsOld_Local[(e-1)*m_QpointsNum+j-1]=m_QpointsScalarMaterials_Local[(e-1)*m_QpointsNum+j-1];
            m_QpointsVectorMaterialsOld_Local[(e-1)*m_QpointsNum+j-1]=m_QpointsVectorMaterials_Local[(e-1)*m_QpointsNum+j-1];
            m_QpointsRank2MaterialsOld_Local[(e-1)*m_QpointsNum+j-1]=m_QpointsRank2Materials_Local[(e-1)*m_QpointsNum+j-1];
            m_QpointsRank4MaterialsOld_Local[(e-1)*m_QpointsNum+j-1]=m_QpointsRank4Materials_Local[(e-1)*m_QpointsNum+j-1];
        }
    }
}