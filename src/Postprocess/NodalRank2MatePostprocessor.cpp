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
//+++ Purpose: Get the nodal rank-2 value via its node id for pps
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/NodalRank2MatePostprocessor.h"

double NodalRank2MatePostprocessor::computeNodalValue(const int &dofid,
                                     const nlohmann::json &parameters,
                                     const DofHandler &dofhandler,
                                     SolutionSystem &soln,
                                     ProjectionSystem &projsystem){
    if(dofid||soln.getBulkElmtsNum()||dofhandler.getBulkElmtsNum()){}

    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"nodeid","rank2mate","i-index","j-index"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the NodalRank2MatePostprocessor, "
                                      "the 'nodeid', 'rank2mate', 'i-index', and 'j-index' are the only parameters you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_nodeid=JsonUtils::getInteger(parameters,"nodeid");
    if(m_nodeid<1||m_nodeid>dofhandler.getNodesNum()){
        MessagePrinter::printErrorTxt("node id="+to_string(m_nodeid)+" is invalid for NodalRank2MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_i=JsonUtils::getInteger(parameters,"i-index");
    if(m_i<1||m_i>3){
        MessagePrinter::printErrorTxt("i-index="+to_string(m_i)+" is invalid for NodalRank2MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_j=JsonUtils::getInteger(parameters,"j-index");
    if(m_j<1||m_j>3){
        MessagePrinter::printErrorTxt("j-index="+to_string(m_j)+" is invalid for NodalRank2MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_rank2matename=JsonUtils::getString(parameters,"rank2mate");
    
    m_pps_value=projsystem.getIthNodeRank2MateViaMateName(m_nodeid,m_rank2matename);
    
    return m_pps_value(m_i,m_j);

}