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
//+++ Date   : 2022.09.30
//+++ Purpose: This class calculates the integration of the specific
//+++          rank-2 material on the specific side
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/SideIntegralRank4MatePostprocessor.h"

double SideIntegralRank4MatePostprocessor::computeSideIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &shp,
                                            SolutionSystem &soln,
                                            ProjectionSystem &projsystem){
    if(dofid||soln.getDofsNum()||elmtinfo.m_Dim){}
    
    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"nodeid","rank4mate","i-index","j-index","k-index","l-index"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the NodalRank4MatePostprocessor, "
                                      "the 'nodeid', 'rank4mate', 'i-index', 'j-index', 'k-index', and 'l-index' are the only parameters you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_i=JsonUtils::getInteger(parameters,"i-index");
    if(m_i<1||m_i>3){
        MessagePrinter::printErrorTxt("i-index="+to_string(m_i)+" is invalid for SideIntegralRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_j=JsonUtils::getInteger(parameters,"j-index");
    if(m_j<1||m_j>3){
        MessagePrinter::printErrorTxt("j-index="+to_string(m_j)+" is invalid for SideIntegralRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_k=JsonUtils::getInteger(parameters,"k-index");
    if(m_k<1||m_k>3){
        MessagePrinter::printErrorTxt("k-index="+to_string(m_k)+" is invalid for SideIntegralRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_l=JsonUtils::getInteger(parameters,"l-index");
    if(m_l<1||m_l>3){
        MessagePrinter::printErrorTxt("l-index="+to_string(m_l)+" is invalid for SideIntegralRank4MatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_rank4matename=JsonUtils::getString(parameters,"rank4mate");
    
    m_rank4value=projsystem.getIthNodeRank4MateViaMateName(nodeid,m_rank4matename);
    
    m_ppsvalue=m_rank4value(m_i,m_j,m_k,m_l)*shp.m_Test;
    
    return m_ppsvalue;
}