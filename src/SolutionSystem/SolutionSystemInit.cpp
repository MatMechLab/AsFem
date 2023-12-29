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
//+++ Purpose: initialize the solution system
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SolutionSystem/SolutionSystem.h"

void SolutionSystem::init(const DofHandler &t_dofhandler,const FE &t_fe){
    m_dofs=t_dofhandler.getActiveDofs();
    m_bulkelmts_num=t_dofhandler.getBulkElmtsNum();
    m_qpoints_num=t_fe.m_bulk_qpoints.getQPointsNum();

    //******************************************************
    //*** initialize each vector
    //******************************************************
    m_u_current.resize(m_dofs,0.0);
    m_u_old.resize(m_dofs,0.0);
    m_u_older.resize(m_dofs,0.0);
    m_u_temp.resize(m_dofs,0.0);
    m_u_copy.resize(m_dofs,0.0);

    // for velocity and acceleration
    m_v.resize(m_dofs,0.0);
    m_a.resize(m_dofs,0.0);

    // for the material properties on each gauss point
    m_qpoints_scalarmaterials.resize(m_bulkelmts_num*m_qpoints_num);
    m_qpoints_vectormaterials.resize(m_bulkelmts_num*m_qpoints_num);
    m_qpoints_rank2materials.resize(m_bulkelmts_num*m_qpoints_num);
    m_qpoints_rank4materials.resize(m_bulkelmts_num*m_qpoints_num);

    m_qpoints_scalarmaterials_old.resize(m_bulkelmts_num*m_qpoints_num);
    m_qpoints_vectormaterials_old.resize(m_bulkelmts_num*m_qpoints_num);
    m_qpoints_rank2materials_old.resize(m_bulkelmts_num*m_qpoints_num);
    m_qpoints_rank4materials_old.resize(m_bulkelmts_num*m_qpoints_num);

    m_allocated=true;

}