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
//+++ Date   : 2022.10.26
//+++ Purpose: This class carry out the side integration for the user-1
//+++          defined calculation.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/User1SideIntegralPostprocessor.h"

double User1SideIntegralPostprocessor::computeSideIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &shp,
                                            SolutionSystem &soln,
                                            ProjectionSystem &projsystem){
    if(nodeid||parameters.size()||projsystem.getNodesNum()||shp.m_test) {}
    if(dofid<1||dofid>soln.getDofsNum()){
        MessagePrinter::printErrorTxt("dof id="+to_string(dofid)+" is out of range for User1SideIntegralPostprocessor");
        MessagePrinter::exitAsFem();
    }
    
    m_ppsvalue=std::cos(elmtinfo.m_gpCoords0(1))/elmtinfo.m_nodesnum;// return cos(x) for int(cos(x))dx integration
    
    return m_ppsvalue;
}