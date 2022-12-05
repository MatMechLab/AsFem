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
//+++ Date   : 2022.07.29
//+++ Purpose: update the solution array and materials array
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SolutionSystem/SolutionSystem.h"

void SolutionSystem::updateSolution(){
    m_u_older=m_u_old;
    m_u_old=m_u_current;
}

void SolutionSystem::updateMaterialsSolution(){
    // TODO: in the future, this vector should be designed as the distributed one (mpi vector) 
    //       to reduce the memory consumption !!!
    for(int e=1;e<=m_bulkelmts_num;e++){
        for(int j=1;j<=m_qpoints_num;j++){
            m_qpoints_scalarmaterials_old[(e-1)*m_qpoints_num+j-1]=m_qpoints_scalarmaterials[(e-1)*m_qpoints_num+j-1];
            m_qpoints_vectormaterials_old[(e-1)*m_qpoints_num+j-1]=m_qpoints_vectormaterials[(e-1)*m_qpoints_num+j-1];
            m_qpoints_rank2materials_old[(e-1)*m_qpoints_num+j-1]=m_qpoints_rank2materials[(e-1)*m_qpoints_num+j-1];
            m_qpoints_rank4materials_old[(e-1)*m_qpoints_num+j-1]=m_qpoints_rank4materials[(e-1)*m_qpoints_num+j-1];
        }
    }
}