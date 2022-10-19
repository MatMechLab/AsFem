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
//+++ Purpose: execute the nodal type postprocess and return a scalar value
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocessor.h"

double Postprocessor::executeNodalPostprocess(const PostprocessorType &pps_type,
                                              const int &dofid,
                                              const nlohmann::json &parameters,
                                              const DofHandler &dofhandler,
                                              SolutionSystem &soln,
                                              ProjectionSystem &projsystem){
    double pps_value=0.0;
    switch (pps_type)
    {
    case PostprocessorType::NODALVALUE:
        pps_value=NodalValuePostprocessor::computeNodalValue(dofid,parameters,dofhandler,soln,projsystem);
        break;
    case PostprocessorType::NODALSCALARMATERIALVALUE:
        pps_value=NodalScalarMatePostprocessor::computeNodalValue(dofid,parameters,dofhandler,soln,projsystem);
        break;
    case PostprocessorType::NODALVECTORMATERIALVALUE:
        pps_value=NodalVectorMatePostprocessor::computeNodalValue(dofid,parameters,dofhandler,soln,projsystem);
        break;
    case PostprocessorType::NODALRANK2MATERIALVALUE:
        pps_value=NodalRank2MatePostprocessor::computeNodalValue(dofid,parameters,dofhandler,soln,projsystem);
        break;
    case PostprocessorType::NODALRANK4MATERIALVALUE:
        pps_value=NodalRank4MatePostprocessor::computeNodalValue(dofid,parameters,dofhandler,soln,projsystem);
        break;
    default:
        MessagePrinter::printErrorTxt("Unsupported nodal value postprocessor in executeNodalPostprocess");
        MessagePrinter::exitAsFem();
        break;
    }
    return pps_value;
}