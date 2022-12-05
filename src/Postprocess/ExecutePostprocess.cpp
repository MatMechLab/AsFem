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
//+++ Date   : 2022.09.23
//+++ Purpose: Execute the postprocess system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocessor.h"

void Postprocessor::executePostprocess(const Mesh &t_mesh,
                                       const DofHandler &t_dofhandler,
                                       FE &t_fe,
                                       MateSystem &t_matesystem,
                                       ProjectionSystem &t_projsystem,
                                       SolutionSystem &t_solution){
    if(t_mesh.getBulkMeshBulkElmtOrder()||t_fe.getMaxDim()||t_matesystem.m_materialcontainer.getRank2MaterialsNum()){}
    
    t_projsystem.makeGhostCopyOfProjectionData();// always make ghost copy here, not in each individual rank !!!
    
    int dofid;
    PostprocessorType pps_type;
    for(int i=1;i<=getPPSBlocksNum();i++){
        dofid=getIthPPSBlock(i).m_dofid;
        pps_type=getIthPPSBlock(i).m_pps_type;
        switch (pps_type)
        {
        case PostprocessorType::NULLPPS:
            break;
        case PostprocessorType::NODALVALUE:
        case PostprocessorType::NODALSCALARMATERIALVALUE:
        case PostprocessorType::NODALVECTORMATERIALVALUE:
        case PostprocessorType::NODALRANK2MATERIALVALUE:
        case PostprocessorType::NODALRANK4MATERIALVALUE:
        {
            m_pps_values[i-1]=executeNodalPostprocess(pps_type,dofid,getIthPPSBlock(i).m_parameters,t_dofhandler,t_solution,t_projsystem);
            break;
        }
        case PostprocessorType::AREA:
        case PostprocessorType::SIDEAVERAGEVALUE:
        case PostprocessorType::SIDEAVERAGESCALARMATERIALVALUE:
        case PostprocessorType::SIDEAVERAGEVECTORMATERIALVALUE:
        case PostprocessorType::SIDEAVERAGERANK2MATERIALVALUE:
        case PostprocessorType::SIDEAVERAGERANK4MATERIALVALUE:
        case PostprocessorType::SIDEINTEGRATEVALUE:
        case PostprocessorType::SIDEINTEGRATESCALARMATERIALVALUE:
        case PostprocessorType::SIDEINTEGRATEVECTORMATERIALVALUE:
        case PostprocessorType::SIDEINTEGRATERANK2MATERIALVALUE:
        case PostprocessorType::SIDEINTEGRATERANK4MATERIALVALUE:
        case PostprocessorType::USER1SIDEINTEGRALPPS:
        {
            m_pps_values[i-1]=executeSideIntegralPostprocess(pps_type,dofid,getIthPPSBlock(i).m_sidenamelist,getIthPPSBlock(i).m_parameters,t_mesh,t_dofhandler,t_fe,t_solution,t_projsystem);
            break;
        }
        case PostprocessorType::VOLUME:
        case PostprocessorType::VOLUMEAVERAGEVALUE:
        case PostprocessorType::VOLUMEAVERAGESCALARMATERIALVALUE:
        case PostprocessorType::VOLUMEAVERAGEVECTORMATERIALVALUE:
        case PostprocessorType::VOLUMEAVERAGERANK2MATERIALVALUE:
        case PostprocessorType::VOLUMEAVERAGERANK4MATERIALVALUE:
        case PostprocessorType::VOLUMEINTEGRATEVALUE:
        case PostprocessorType::VOLUMEINTEGRATESCALARMATERIALVALUE:
        case PostprocessorType::VOLUMEINTEGRATEVECTORMATERIALVALUE:
        case PostprocessorType::VOLUMEINTEGRATERANK2MATERIALVALUE:
        case PostprocessorType::VOLUMEINTEGRATERANK4MATERIALVALUE:
        case PostprocessorType::USER1VOLUMEINTEGRALPPS:
        {
            m_pps_values[i-1]=executeVolumeIntegralPostprocess(pps_type,dofid,getIthPPSBlock(i).m_domainnamelist,getIthPPSBlock(i).m_parameters,t_mesh,t_dofhandler,t_fe,t_solution,t_projsystem);
            break;
        }
        default:
            MessagePrinter::printErrorTxt("Unsupported postprocessor type in executePostprocess, please check either your code or your input file");
            MessagePrinter::exitAsFem();
            break;
        }
    }
    t_projsystem.destroyGhostCopyOfProjectionData();
}