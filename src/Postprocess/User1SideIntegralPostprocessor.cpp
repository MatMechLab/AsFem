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
    if(nodeid||parameters.size()||projsystem.getNodesNum()||shp.m_Test) {}
    if(dofid<1||dofid>soln.getDofsNum()){
        MessagePrinter::printErrorTxt("dof id="+to_string(dofid)+" is out of range for User1SideIntegralPostprocessor");
        MessagePrinter::exitAsFem();
    }
    
    m_ppsvalue=std::cos(elmtinfo.m_QpCoords0(1))/elmtinfo.m_NodesNum;// return cos(x) for int(cos(x))dx integration
    
    return m_ppsvalue;
}