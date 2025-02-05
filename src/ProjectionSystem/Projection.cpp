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
//+++ Date   : 2022.08.22
//+++ Purpose: execute the projection from integration points to
//+++          nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ProjectionSystem/ProjectionSystem.h"

void ProjectionSystem::executeProjection(const FECell &t_FECell,
                                         const DofHandler &t_DofHandler,
                                         const ElmtSystem &t_ElmtSystem,
                                         MateSystem &t_MateSystem,
                                         FE &t_FE,
                                         SolutionSystem &t_SolnSystem,
                                         const FEControlInfo &t_FECtrlInfo){
    // if no projection is required, then return back
    if(getScalarMaterialNum()<1 &&
       getVectorMaterialNum()<1 &&
       getRank2MaterialNum()<1 &&
       getRank4MaterialNum()<1){
        return;
    }
    if (m_ProjType==ProjectionType::DEFAULT||
        m_ProjType==ProjectionType::LEASTSQUARE) {
        LeastSquareProjection::executeMyProjection(t_FECell,t_DofHandler,t_ElmtSystem,t_MateSystem,t_FE,t_SolnSystem,t_FECtrlInfo,m_Data);
    }
    else if (m_ProjType==ProjectionType::FULLLEASTSQUARE) {
        FullLeastSquareProjection::executeMyProjection(t_FECell,t_DofHandler,t_ElmtSystem,t_MateSystem,t_FE,t_SolnSystem,t_FECtrlInfo,m_Data);
    }
}