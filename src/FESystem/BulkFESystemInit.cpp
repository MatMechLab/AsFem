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
//+++ Date   : 2022.07.29
//+++ Purpose: initialize the bulk fe system
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/BulkFESystem.h"

void BulkFESystem::init(const FECell &t_fecell,const DofHandler &t_dofhandler){
    m_MaxElmtDofs=t_dofhandler.getMaxDofsPerElmt();
    m_MaxNodalDofs=t_dofhandler.getMaxDofsPerNode();
    m_BulkElmtNodesNum=t_fecell.getFECellNodesNumPerBulkElmt();

    m_LocalElmtInfo.m_Dim=t_fecell.getFECellMaxDim();
    m_LocalElmtInfo.m_NodesNum=m_BulkElmtNodesNum;

    m_MaxKmatCoeff=-1.0e16;

    m_LocalR.resize(m_MaxElmtDofs+1,0.0);
    m_LocalK.resize(m_MaxElmtDofs+1,m_MaxElmtDofs+1,0.0);

    m_SubR.resize(m_MaxNodalDofs+1,0.0);
    m_SubK.resize(m_MaxNodalDofs+1,m_MaxNodalDofs+1,0.0);

    m_ElmtConn.resize(m_BulkElmtNodesNum,0);
    m_ElmtDofIDs.resize(m_MaxElmtDofs+1,0);
    m_SubElmtDofIDs.resize(m_MaxElmtDofs+1,0);

    // for elemental solution
    m_ElmtU.resize(m_MaxElmtDofs,0.0);
    m_ElmtUold.resize(m_MaxElmtDofs,0.0);
    m_ElmtUolder.resize(m_MaxElmtDofs,0.0);

    m_ElmtV.resize(m_MaxElmtDofs,0.0);
    m_ElmtA.resize(m_MaxElmtDofs,0.0);

    // for sub elemental solution
    m_LocalElmtSoln.m_QpU.resize(m_MaxNodalDofs+1,0.0);
    m_LocalElmtSoln.m_QpUold.resize(m_MaxNodalDofs+1,0.0);
    m_LocalElmtSoln.m_QpUolder.resize(m_MaxNodalDofs+1,0.0);

    m_LocalElmtSoln.m_QpV.resize(m_MaxNodalDofs+1,0.0);
    m_LocalElmtSoln.m_QpA.resize(m_MaxNodalDofs+1,0.0);

    m_LocalElmtSoln.m_QpGradU.resize(m_MaxNodalDofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_QpGradUold.resize(m_MaxNodalDofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_QpGradUolder.resize(m_MaxNodalDofs+1,Vector3d(0.0));

    m_LocalElmtSoln.m_QpGradV.resize(m_MaxNodalDofs+1,Vector3d(0.0));

    m_LocalElmtSoln.m_Qpgradu.resize(m_MaxNodalDofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_Qpgradv.resize(m_MaxNodalDofs+1,Vector3d(0.0));


    m_Nodes.resize(m_BulkElmtNodesNum);
    m_Nodes0.resize(m_BulkElmtNodesNum);

}