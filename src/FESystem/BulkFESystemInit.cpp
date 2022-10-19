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
//+++ Date   : 2022.07.29
//+++ Purpose: initialize the bulk fe system
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/BulkFESystem.h"

void BulkFESystem::init(const Mesh &t_mesh,const DofHandler &t_dofhandler){
    m_max_elmt_dofs=t_dofhandler.getMaxDofsPerElmt();
    m_max_nodal_dofs=t_dofhandler.getMaxDofsPerNode();
    m_bulkelmt_nodesnum=t_mesh.getBulkMeshNodesNumPerBulkElmt();

    m_local_elmtinfo.m_dim=t_mesh.getBulkMeshMaxDim();
    m_local_elmtinfo.m_nodesnum=m_bulkelmt_nodesnum;

    m_max_k_coeff=-1.0e16;

    m_localR.resize(m_max_elmt_dofs+1,0.0);
    m_localK.resize(m_max_elmt_dofs+1,m_max_elmt_dofs+1,0.0);

    m_subR.resize(m_max_nodal_dofs+1,0.0);
    m_subK.resize(m_max_nodal_dofs+1,m_max_nodal_dofs+1,0.0);

    m_elmtconn.resize(m_bulkelmt_nodesnum,0);
    m_elmtdofsid.resize(m_max_elmt_dofs+1,0);
    m_subelmtdofsid.resize(m_max_nodal_dofs+1,0);

    // for elemental solution
    m_elmtU.resize(m_max_elmt_dofs,0.0);
    m_elmtUold.resize(m_max_elmt_dofs,0.0);
    m_elmtUolder.resize(m_max_elmt_dofs,0.0);

    m_elmtV.resize(m_max_elmt_dofs,0.0);
    m_elmtA.resize(m_max_elmt_dofs,0.0);

    // for sub elemental solution
    m_local_elmtsoln.m_gpU.resize(m_max_nodal_dofs+1,0.0);
    m_local_elmtsoln.m_gpUold.resize(m_max_nodal_dofs+1,0.0);
    m_local_elmtsoln.m_gpUolder.resize(m_max_nodal_dofs+1,0.0);

    m_local_elmtsoln.m_gpV.resize(m_max_nodal_dofs+1,0.0);
    m_local_elmtsoln.m_gpA.resize(m_max_nodal_dofs+1,0.0);

    m_local_elmtsoln.m_gpGradU.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_local_elmtsoln.m_gpGradUold.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_local_elmtsoln.m_gpGradUolder.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_local_elmtsoln.m_gpGradV.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_local_elmtsoln.m_gpgradu.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_local_elmtsoln.m_gpgradv.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_local_elmtinfo.m_gpCoords=0.0;
    m_local_elmtinfo.m_gpCoords0=0.0;

    m_nodes.resize(m_bulkelmt_nodesnum);
    m_nodes0.resize(m_bulkelmt_nodesnum);

}