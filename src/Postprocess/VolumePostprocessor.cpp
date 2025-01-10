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
//+++ Date   : 2022.09.29
//+++ Purpose: This class calculates the area of the specific side
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/VolumePostprocessor.h"

double VolumePostprocessor::computeVolumeIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &t_parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &t_shp,
                                            SolutionSystem &t_soln,
                                            ProjectionSystem &t_projsystem){
    if(dofid||nodeid||t_parameters.size()||t_shp.m_Test||t_soln.getDofsNum()||t_projsystem.getNodesNum()) {}
    
    
    m_ppsvalue=1.0/elmtinfo.m_NodesNum;
    
    return m_ppsvalue;
}