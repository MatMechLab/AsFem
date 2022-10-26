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
//+++          scalar material on the specific side
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/SideIntegralScalarMatePostprocessor.h"

double SideIntegralScalarMatePostprocessor::computeSideIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &shp,
                                            SolutionSystem &soln,
                                            ProjectionSystem &projsystem){
    if(dofid||soln.getDofsNum()||elmtinfo.m_dim){}
    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"scalarmate"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the SideIntegralScalarMatePostprocessor, "
                                      "the 'scalarmate' is the only parameter you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_scalarmatename=JsonUtils::getString(parameters,"scalarmate");
    
    m_scalarvalue=projsystem.getIthNodeScalarMateViaMateName(nodeid,m_scalarmatename);
    
    m_ppsvalue=m_scalarvalue*shp.m_test;
    
    return m_ppsvalue;
}