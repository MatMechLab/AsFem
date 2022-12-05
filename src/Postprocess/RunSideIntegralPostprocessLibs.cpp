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
//+++ Purpose: execute the side integral postprocess and return a scalar value
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocessor.h"

double Postprocessor::runSideIntegralPostprocessLibs(const PostprocessorType &pps_type,
                                          const int &dofid,
                                          const int &nodeid,
                                          const nlohmann::json &t_parameters,
                                          const LocalElmtInfo &t_elmtinfo,
                                          const LocalShapeFun &t_shp,
                                          SolutionSystem &t_soln,
                                          ProjectionSystem &t_projsystem){
    double pps_result;
    pps_result=0.0;
    switch (pps_type)
    {
    case PostprocessorType::AREA:
        pps_result=AreaPostprocessor::computeSideIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::SIDEINTEGRATEVALUE:
    case PostprocessorType::SIDEAVERAGEVALUE:
        pps_result=SideIntegralValuePostprocessor::computeSideIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::SIDEINTEGRATESCALARMATERIALVALUE:
    case PostprocessorType::SIDEAVERAGESCALARMATERIALVALUE:
        pps_result=SideIntegralScalarMatePostprocessor::computeSideIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::SIDEINTEGRATEVECTORMATERIALVALUE:
    case PostprocessorType::SIDEAVERAGEVECTORMATERIALVALUE:
        pps_result=SideIntegralVectorMatePostprocessor::computeSideIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::SIDEINTEGRATERANK2MATERIALVALUE:
    case PostprocessorType::SIDEAVERAGERANK2MATERIALVALUE:
        pps_result=SideIntegralRank2MatePostprocessor::computeSideIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::SIDEINTEGRATERANK4MATERIALVALUE:
    case PostprocessorType::SIDEAVERAGERANK4MATERIALVALUE:
        pps_result=SideIntegralRank4MatePostprocessor::computeSideIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::USER1SIDEINTEGRALPPS:
        pps_result=User1SideIntegralPostprocessor::computeSideIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    default:
        MessagePrinter::printErrorTxt("Unsupported side integral postprocess type in runSideIntegralPostprocessLibs, "
                                      "please check your code or your input file");
        break;
        
    }
    return pps_result;
}