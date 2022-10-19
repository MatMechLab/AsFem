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
//+++ Date   : 2022.06.13
//+++ Purpose: defines the solution array for FEM analysis, stores
//+++          the necessary results.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SolutionSystem/SolutionSystem.h"

SolutionSystem::SolutionSystem(){
    m_dofs=0;
    m_allocated=false;
    m_bulkelmts_num=0;
    m_qpoints_num=0;
}

void SolutionSystem::releaseMemory(){
    if(m_allocated){
        m_u_current.releaseMemory();
        m_u_old.releaseMemory();
        m_u_older.releaseMemory();
        m_u_temp.releaseMemory();
        m_u_copy.releaseMemory();

        m_v.releaseMemory();

        m_a.releaseMemory();

        for(auto &it:m_qpoints_scalarmaterials) it.clear();
        m_qpoints_scalarmaterials.clear();

        for(auto &it:m_qpoints_vectormaterials) it.clear();
        m_qpoints_vectormaterials.clear();

        for(auto &it:m_qpoints_rank2materials) it.clear();
        m_qpoints_rank2materials.clear();

        for(auto &it:m_qpoints_rank4materials) it.clear();
        m_qpoints_rank4materials.clear();

        //****************************************
        for(auto &it:m_qpoints_scalarmaterials_old) it.clear();
        m_qpoints_scalarmaterials_old.clear();

        for(auto &it:m_qpoints_vectormaterials_old) it.clear();
        m_qpoints_vectormaterials_old.clear();

        for(auto &it:m_qpoints_rank2materials_old) it.clear();
        m_qpoints_rank2materials_old.clear();

        for(auto &it:m_qpoints_rank4materials_old) it.clear();
        m_qpoints_rank4materials_old.clear();

        m_allocated=false;
    }
}