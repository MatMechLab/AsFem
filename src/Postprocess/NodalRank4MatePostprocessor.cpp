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
//+++ Purpose: Get the nodal rank-4 material value via its node id for pps
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/NodalRank4MatePostprocessor.h"

double NodalRank4MatePostprocessor::computeNodalValue(const int &dofid,
                                     const nlohmann::json &parameters,
                                     const DofHandler &dofhandler,
                                     SolutionSystem &soln,
                                     ProjectionSystem &projsystem){
    if(dofid||soln.getBulkElmtsNum()||dofhandler.getBulkElmtsNum()){}

    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"nodeid","rank4mate","i-index","j-index","k-index","l-index"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the NodalRank4MatePostprocessor, "
                                      "the 'nodeid', 'rank4mate', 'i-index', 'j-index', 'k-index', and 'l-index' are the only parameters you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_nodeid=JsonUtils::getInteger(parameters,"nodeid");
    if(m_nodeid<1||m_nodeid>dofhandler.getNodesNum()){
        MessagePrinter::printErrorTxt("node id="+to_string(m_nodeid)+" is invalid for NodalRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_i=JsonUtils::getInteger(parameters,"i-index");
    if(m_i<1||m_i>3){
        MessagePrinter::printErrorTxt("i-index="+to_string(m_i)+" is invalid for NodalRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_j=JsonUtils::getInteger(parameters,"j-index");
    if(m_j<1||m_j>3){
        MessagePrinter::printErrorTxt("j-index="+to_string(m_j)+" is invalid for NodalRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_k=JsonUtils::getInteger(parameters,"k-index");
    if(m_k<1||m_k>3){
        MessagePrinter::printErrorTxt("k-index="+to_string(m_k)+" is invalid for NodalRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_l=JsonUtils::getInteger(parameters,"l-index");
    if(m_l<1||m_l>3){
        MessagePrinter::printErrorTxt("l-index="+to_string(m_l)+" is invalid for NodalRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_rank4matename=JsonUtils::getString(parameters,"rank4mate");
    
    m_pps_value=projsystem.getIthNodeRank4MateViaMateName(m_nodeid,m_rank4matename);
    
    return m_pps_value(m_i,m_j,m_k,m_l);

}