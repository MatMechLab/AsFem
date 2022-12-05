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
//+++          vector material on the specific side
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/SideIntegralVectorMatePostprocessor.h"

double SideIntegralVectorMatePostprocessor::computeSideIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &shp,
                                            SolutionSystem &soln,
                                            ProjectionSystem &projsystem){
    if(dofid||soln.getDofsNum()||elmtinfo.m_dim){}
    
    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"vectormate","component"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the SideIntegralVectorMatePostprocessor, "
                                      "the 'vectormate', and 'component' are the only parameters you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_component=JsonUtils::getInteger(parameters,"component");
    if(m_component<1||m_component>3){
        MessagePrinter::printErrorTxt("component="+to_string(m_component)+" is invalid for SideIntegralVectorMatePostprocessor");
        MessagePrinter::exitAsFem();
    }

    m_vectormatename=JsonUtils::getString(parameters,"vectormate");
    
    m_vectorvalue=projsystem.getIthNodeVectorMateViaMateName(nodeid,m_vectormatename);
    
    m_ppsvalue=m_vectorvalue(m_component)*shp.m_test;
    
    return m_ppsvalue;
}