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
//+++ Purpose: Get the nodal scalar material value via its node id for pps
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/NodalScalarMatePostprocessor.h"

double NodalScalarMatePostprocessor::computeNodalValue(const int &dofid,
                                     const nlohmann::json &parameters,
                                     const DofHandler &dofhandler,
                                     SolutionSystem &soln,
                                     ProjectionSystem &projsystem){
    if(dofid||soln.getDofsNum()||dofhandler.getBulkElmtsNum()){}

    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"nodeid","scalarmate"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the NodalScalarMatePostprocessor, "
                                      "the 'nodeid' and 'scalarmate' are the only parameters you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_nodeid=JsonUtils::getInteger(parameters,"nodeid");
    if(m_nodeid<1||m_nodeid>dofhandler.getNodesNum()){
        MessagePrinter::printErrorTxt("node id="+to_string(m_nodeid)+" is invalid for NodalScalarMatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_scalarmatename=JsonUtils::getString(parameters,"scalarmate");
    
    m_pps_value=projsystem.getIthNodeScalarMateViaMateName(m_nodeid,m_scalarmatename);
    
    return m_pps_value;

}