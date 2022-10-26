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
//+++ Date   : 2022.09.30
//+++ Purpose: This class calculates the integration of the specific
//+++          rank-2 material on the specific side
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/SideIntegralRank2MatePostprocessor.h"

double SideIntegralRank2MatePostprocessor::computeSideIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &shp,
                                            SolutionSystem &soln,
                                            ProjectionSystem &projsystem){
    if(dofid||soln.getDofsNum()||elmtinfo.m_dim){}
    
    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"rank2mate","i-index","j-index"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the SideIntegralRank2MatePostprocessor, "
                                      "the 'rank2mate', 'i-index', and 'j-index' are the only parameters you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_i=JsonUtils::getInteger(parameters,"i-index");
    if(m_i<1||m_i>3){
        MessagePrinter::printErrorTxt("i-index="+to_string(m_i)+" is invalid for SideIntegralRank2MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_j=JsonUtils::getInteger(parameters,"j-index");
    if(m_j<1||m_j>3){
        MessagePrinter::printErrorTxt("j-index="+to_string(m_j)+" is invalid for SideIntegralRank2MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_rank2matename=JsonUtils::getString(parameters,"rank2mate");
    
    m_rank2value=projsystem.getIthNodeRank2MateViaMateName(nodeid,m_rank2matename);
    
    m_ppsvalue=m_rank2value(m_i,m_j)*shp.m_test;
    
    return m_ppsvalue;
}