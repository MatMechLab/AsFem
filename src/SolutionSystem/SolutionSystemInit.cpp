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
    m_Dofs=t_dofhandler.getActiveDofs();
    m_BulkElmtsNum=t_dofhandler.getBulkElmtsNum();
    m_BulkElmtsNum_Local=t_dofhandler.getLocalBulkElmtsNum();
    m_QpointsNum=t_fe.m_BulkQpoints.getQPointsNum();

    //******************************************************
    //*** initialize each vector
    //******************************************************
    m_Ucurrent.resize(m_Dofs,0.0);
    m_Uold.resize(m_Dofs,0.0);
    m_Uolder.resize(m_Dofs,0.0);
    m_Utemp.resize(m_Dofs,0.0);
    m_Ucopy.resize(m_Dofs,0.0);
    m_dU.resize(m_Dofs,0.0);

    // for velocity and acceleration
    m_V.resize(m_Dofs,0.0);
    m_A.resize(m_Dofs,0.0);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    

    if(rank==0){
        // for the material properties on each gauss point
        m_QpointsScalarMaterials_Total.resize(m_BulkElmtsNum*m_QpointsNum);
        m_QpointsVectorMaterials_Total.resize(m_BulkElmtsNum*m_QpointsNum);
        m_QpointsRank2Materials_Total.resize(m_BulkElmtsNum*m_QpointsNum);
        m_QpointsRank4Materials_Total.resize(m_BulkElmtsNum*m_QpointsNum);
        // for the old materials
        m_QpointsScalarMaterialsOld_Total.resize(m_BulkElmtsNum*m_QpointsNum);
        m_QpointsVectorMaterialsOld_Total.resize(m_BulkElmtsNum*m_QpointsNum);
        m_QpointsRank2MaterialsOld_Total.resize(m_BulkElmtsNum*m_QpointsNum);
        m_QpointsRank4MaterialsOld_Total.resize(m_BulkElmtsNum*m_QpointsNum);
    }

    // for the material properties on each gauss point owned by each rank
    m_QpointsScalarMaterials_Local.resize(m_BulkElmtsNum_Local*m_QpointsNum);
    m_QpointsVectorMaterials_Local.resize(m_BulkElmtsNum_Local*m_QpointsNum);
    m_QpointsRank2Materials_Local.resize(m_BulkElmtsNum_Local*m_QpointsNum);
    m_QpointsRank4Materials_Local.resize(m_BulkElmtsNum_Local*m_QpointsNum);
    // for the old materials
    m_QpointsScalarMaterialsOld_Local.resize(m_BulkElmtsNum_Local*m_QpointsNum);
    m_QpointsVectorMaterialsOld_Local.resize(m_BulkElmtsNum_Local*m_QpointsNum);
    m_QpointsRank2MaterialsOld_Local.resize(m_BulkElmtsNum_Local*m_QpointsNum);
    m_QpointsRank4MaterialsOld_Local.resize(m_BulkElmtsNum_Local*m_QpointsNum);

    m_Allocated=true;

}