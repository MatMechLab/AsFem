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
//+++ Purpose: Get the nodal vector material value via its node id for pps
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/NodalVectorMatePostprocessor.h"

double NodalVectorMatePostprocessor::computeNodalValue(const int &dofid,
                                     const nlohmann::json &parameters,
                                     const DofHandler &dofhandler,
                                     SolutionSystem &soln,
                                     ProjectionSystem &projsystem){
    if(dofid||soln.getDofsNum()||dofhandler.getBulkElmtsNum()){}

    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"nodeid","vectormate","component"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the NodalVectorMatePostprocessor, "
                                      "the 'nodeid', 'vectormate', and 'component' are the only parameters you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_nodeid=JsonUtils::getInteger(parameters,"nodeid");
    if(m_nodeid<1||m_nodeid>dofhandler.getNodesNum()){
        MessagePrinter::printErrorTxt("node id="+to_string(m_nodeid)+" is invalid for NodalVectorMatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_component=JsonUtils::getInteger(parameters,"component");
    if(m_component<1||m_component>3){
        MessagePrinter::printErrorTxt("component="+to_string(m_component)+" is invalid for NodalVectorMatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_vectormatename=JsonUtils::getString(parameters,"vectormate");
    
    m_pps_value=projsystem.getIthNodeVectorMateViaMateName(m_nodeid,m_vectormatename);
    
    return m_pps_value(m_component);

}