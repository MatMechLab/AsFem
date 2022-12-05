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
//+++ Date   : 2022.09.29
//+++ Purpose: This class calculates the area of the specific side
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/AreaPostprocessor.h"

double AreaPostprocessor::computeSideIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &shp,
                                            SolutionSystem &soln,
                                            ProjectionSystem &projsystem){
    if(dofid||nodeid||parameters.size()||shp.m_test||soln.getDofsNum()||projsystem.getNodesNum()) {}
    
    m_ppsvalue=1.0/elmtinfo.m_nodesnum;
    
    return m_ppsvalue;
}