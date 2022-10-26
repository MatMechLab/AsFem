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
//+++ Date   : 2022.09.28
//+++ Purpose: This class carry out the side integral for the specific
//+++          dof, and return the final ingetrated result(scalar)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/SideIntegralValuePostprocessor.h"

double SideIntegralValuePostprocessor::computeSideIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &shp,
                                            SolutionSystem &soln,
                                            ProjectionSystem &projsystem){
    if(nodeid||parameters.size()||projsystem.getNodesNum()||elmtinfo.m_dim) {}
    if(dofid<1||dofid>soln.getDofsNum()){
        MessagePrinter::printErrorTxt("dof id="+to_string(dofid)+" is out of range for SideIntegralValuePostprocessor");
        MessagePrinter::exitAsFem();
    }
    
    m_ppsvalue=soln.m_u_current.getIthValueFromGhost(dofid)*shp.m_test;
    
    return m_ppsvalue;
}