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
//+++ Date   : 2022.09.28
//+++ Purpose: Get the nodal value via its node id for pps
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/NodalValuePostprocessor.h"

double NodalValuePostprocessor::computeNodalValue(const int &dofid,
                                     const nlohmann::json &parameters,
                                     const DofHandler &dofhandler,
                                     SolutionSystem &soln,
                                     ProjectionSystem &projsystem){
    if(projsystem.getNodesNum()){}
    if(dofid<1||dofid>dofhandler.getMaxDofsPerNode()){
        MessagePrinter::printErrorTxt("dof id="+to_string(dofid)+" is invalid for NodalValuePostprocessor");
        MessagePrinter::exitAsFem();
    }
    
    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"nodeid"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the NodalValuePostprocessor, "
                                      "the 'nodeid' is the only parameter you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_nodeid=JsonUtils::getInteger(parameters,"nodeid");
    int iInd=dofhandler.getIthNodeJthDofID(m_nodeid,dofid);
    soln.m_Ucurrent.makeGhostCopy();
    m_pps_value=soln.m_Ucurrent.getIthValueFromGhost(iInd);
    soln.m_Ucurrent.destroyGhostCopy();
    return m_pps_value;

}