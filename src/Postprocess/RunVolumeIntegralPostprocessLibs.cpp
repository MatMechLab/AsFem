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
//+++ Date   : 2022.10.10
//+++ Purpose: execute the volume integral postprocess and return a scalar value
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocessor.h"

double Postprocessor::runVolumeIntegralPostprocessLibs(const PostprocessorType &pps_type,
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
    case PostprocessorType::VOLUME:
        pps_result=VolumePostprocessor::computeVolumeIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::VOLUMEINTEGRATEVALUE:
    case PostprocessorType::VOLUMEAVERAGEVALUE:
        pps_result=VolumeIntegralValuePostprocessor::computeVolumeIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::VOLUMEINTEGRATESCALARMATERIALVALUE:
    case PostprocessorType::VOLUMEAVERAGESCALARMATERIALVALUE:
        pps_result=VolumeIntegralScalarMatePostprocessor::computeVolumeIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::VOLUMEINTEGRATEVECTORMATERIALVALUE:
    case PostprocessorType::VOLUMEAVERAGEVECTORMATERIALVALUE:
        pps_result=VolumeIntegralVectorMatePostprocessor::computeVolumeIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::VOLUMEINTEGRATERANK2MATERIALVALUE:
    case PostprocessorType::VOLUMEAVERAGERANK2MATERIALVALUE:
        pps_result=VolumeIntegralRank2MatePostprocessor::computeVolumeIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::VOLUMEINTEGRATERANK4MATERIALVALUE:
    case PostprocessorType::VOLUMEAVERAGERANK4MATERIALVALUE:
        pps_result=VolumeIntegralRank4MatePostprocessor::computeVolumeIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    case PostprocessorType::USER1VOLUMEINTEGRALPPS:
        pps_result=User1VolumeIntegralPostprocessor::computeVolumeIntegralValue(dofid,nodeid,t_parameters,t_elmtinfo,t_shp,t_soln,t_projsystem);
        break;
    default:
        MessagePrinter::printErrorTxt("Unsupported side integral postprocess type in runVolumeIntegralPostprocessLibs, "
                                      "please check your code or your input file");
        break;
        
    }
    return pps_result;
}