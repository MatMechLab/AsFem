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
//+++ Date   : 2022.10.10
//+++ Purpose: This class calculates the volume integral for the
//+++          specific scalar material on specific domain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/VolumeIntegralScalarMatePostprocessor.h"

double VolumeIntegralScalarMatePostprocessor::computeVolumeIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &shp,
                                            SolutionSystem &soln,
                                            ProjectionSystem &projsystem){
    if(dofid||soln.getDofsNum()||projsystem.getNodesNum()||elmtinfo.m_Dim) {}
    
    if(!JsonUtils::hasOnlyGivenValues(parameters,vector<string>{"scalarmate"})){
        MessagePrinter::printErrorTxt("Unsupported options in parameters of the VolumeIntegralScalarMatePostprocessor, "
                                      "the 'scalarmate' is the only parameter you need,"
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_scalarmatename=JsonUtils::getString(parameters,"scalarmate");
    
    m_scalarvalue=projsystem.getIthNodeScalarMateViaMateName(nodeid,m_scalarmatename);
    
    m_ppsvalue=m_scalarvalue*shp.m_Test;
    
    return m_ppsvalue;
}