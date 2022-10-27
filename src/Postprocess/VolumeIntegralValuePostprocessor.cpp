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
//+++ Purpose: This class calculates the volume integral for the
//+++          specific dof on specific domain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/VolumeIntegralValuePostprocessor.h"

double VolumeIntegralValuePostprocessor::computeVolumeIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &t_parameters,
                                            const LocalElmtInfo &elmtinfo,
                                            const LocalShapeFun &t_shp,
                                            SolutionSystem &t_soln,
                                            ProjectionSystem &t_projsystem){
    if(nodeid||t_parameters.size()||t_projsystem.getNodesNum()||elmtinfo.m_dim) {}
    
    if(dofid<1||dofid>t_soln.getDofsNum()){
        MessagePrinter::printErrorTxt("dof id="+to_string(dofid)+" is out of range for VolumeIntegralValuePostprocessor");
        MessagePrinter::exitAsFem();
    }
    
    m_ppsvalue=t_soln.m_u_current.getIthValueFromGhost(dofid)*t_shp.m_test;
    
    return m_ppsvalue;
}