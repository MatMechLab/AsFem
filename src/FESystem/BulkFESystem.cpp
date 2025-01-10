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
//+++ Date   : 2020.11.29
//+++ Purpose: Implement the general tasks of FEM calculation in AsFem,
//+++          i.e. compute residual, compute jacobian
//+++          projection from gauss point to nodal point
//+++          assemble from local element to global, ...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/BulkFESystem.h"

BulkFESystem::BulkFESystem(){
    m_MaxNodalDofs=0;
    m_MaxElmtDofs=0;
    m_BulkElmtNodesNum=0;

    m_LocalElmtInfo.m_Dim=0;
    m_LocalElmtInfo.m_DofsNum=0;
    m_LocalElmtInfo.m_NodesNum=0;
    m_LocalElmtInfo.m_T=0.0;
    m_LocalElmtInfo.m_Dt=0.0;

    m_LocalR.clean();
    m_LocalK.clean();
    m_SubR.clean();
    m_SubK.clean();

    m_ElmtConn.clear();
    m_ElmtDofIDs.clear();
    m_SubElmtDofIDs.clear();

    m_MaxKmatCoeff=-1.0e16;

    m_ElmtU.clear();
    m_ElmtUold.clear();
    m_ElmtUolder.clear();
    m_ElmtV.clear();
    m_ElmtA.clear();

    m_LocalElmtSoln.m_QpU.clear();
    m_LocalElmtSoln.m_QpUold.clear();
    m_LocalElmtSoln.m_QpUolder.clear();
    m_LocalElmtSoln.m_QpV.clear();
    m_LocalElmtSoln.m_QpA.clear();

    m_LocalElmtSoln.m_QpGradU.clear();
    m_LocalElmtSoln.m_QpGradUold.clear();
    m_LocalElmtSoln.m_QpGradUolder.clear();

    m_LocalElmtSoln.m_QpGradV.clear();

    m_LocalElmtSoln.m_Qpgradu.clear();
    m_LocalElmtSoln.m_Qpgradv.clear();

    m_LocalElmtInfo.m_QpCoords0=0.0;
    m_LocalElmtInfo.m_QpCoords =0.0;

    m_Nodes.clear();
    m_Nodes0.clear();

}

void BulkFESystem::releaseMemory(){
    m_MaxNodalDofs=0;
    m_MaxElmtDofs=0;
    m_BulkElmtNodesNum=0;

    m_LocalElmtInfo.m_Dim=0;
    m_LocalElmtInfo.m_DofsNum=0;
    m_LocalElmtInfo.m_NodesNum=0;
    m_LocalElmtInfo.m_T=0.0;
    m_LocalElmtInfo.m_Dt=0.0;

    m_LocalR.clean();
    m_LocalK.clean();
    m_SubR.clean();
    m_SubK.clean();

    m_ElmtConn.clear();
    m_ElmtDofIDs.clear();
    m_SubElmtDofIDs.clear();

    m_MaxKmatCoeff=-1.0e16;

    m_ElmtU.clear();
    m_ElmtUold.clear();
    m_ElmtUolder.clear();
    m_ElmtV.clear();
    m_ElmtA.clear();

    m_LocalElmtSoln.m_QpU.clear();
    m_LocalElmtSoln.m_QpUold.clear();
    m_LocalElmtSoln.m_QpUolder.clear();
    m_LocalElmtSoln.m_QpV.clear();
    m_LocalElmtSoln.m_QpA.clear();

    m_LocalElmtSoln.m_QpGradU.clear();
    m_LocalElmtSoln.m_QpGradUold.clear();
    m_LocalElmtSoln.m_QpGradUolder.clear();

    m_LocalElmtSoln.m_QpGradV.clear();

    m_LocalElmtSoln.m_Qpgradu.clear();
    m_LocalElmtSoln.m_Qpgradv.clear();

    m_Nodes.clear();
    m_Nodes0.clear();
}