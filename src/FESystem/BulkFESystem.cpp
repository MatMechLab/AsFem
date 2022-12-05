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
//+++ Date   : 2020.11.29
//+++ Purpose: Implement the general tasks of FEM calculation in AsFem,
//+++          i.e. compute residual, compute jacobian
//+++          projection from gauss point to nodal point
//+++          assemble from local element to global, ...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/BulkFESystem.h"

BulkFESystem::BulkFESystem(){
    m_max_nodal_dofs=0;
    m_max_elmt_dofs=0;
    m_bulkelmt_nodesnum=0;

    m_local_elmtinfo.m_dim=0;
    m_local_elmtinfo.m_dofsnum=0;
    m_local_elmtinfo.m_nodesnum=0;
    m_local_elmtinfo.m_t=0.0;
    m_local_elmtinfo.m_dt=0.0;

    m_localR.clean();
    m_localK.clean();
    m_subR.clean();
    m_subK.clean();

    m_elmtconn.clear();
    m_elmtdofsid.clear();
    m_subelmtdofsid.clear();

    m_max_k_coeff=-1.0e16;

    m_elmtU.clear();
    m_elmtUold.clear();
    m_elmtUolder.clear();
    m_elmtV.clear();
    m_elmtA.clear();

    m_local_elmtsoln.m_gpU.clear();
    m_local_elmtsoln.m_gpUold.clear();
    m_local_elmtsoln.m_gpUolder.clear();
    m_local_elmtsoln.m_gpV.clear();
    m_local_elmtsoln.m_gpA.clear();

    m_local_elmtsoln.m_gpGradU.clear();
    m_local_elmtsoln.m_gpGradUold.clear();
    m_local_elmtsoln.m_gpGradUolder.clear();

    m_local_elmtsoln.m_gpGradV.clear();

    m_local_elmtsoln.m_gpgradu.clear();
    m_local_elmtsoln.m_gpgradv.clear();

    m_local_elmtinfo.m_gpCoords0=0.0;
    m_local_elmtinfo.m_gpCoords =0.0;

    m_nodes.clear();
    m_nodes0.clear();

}

void BulkFESystem::releaseMemory(){
    m_max_nodal_dofs=0;
    m_max_elmt_dofs=0;
    m_bulkelmt_nodesnum=0;

    m_local_elmtinfo.m_dim=0;
    m_local_elmtinfo.m_dofsnum=0;
    m_local_elmtinfo.m_nodesnum=0;
    m_local_elmtinfo.m_t=0.0;
    m_local_elmtinfo.m_dt=0.0;

    m_localR.clean();
    m_localK.clean();
    m_subR.clean();
    m_subK.clean();

    m_elmtconn.clear();
    m_elmtdofsid.clear();
    m_subelmtdofsid.clear();

    m_max_k_coeff=-1.0e16;

    m_elmtU.clear();
    m_elmtUold.clear();
    m_elmtUolder.clear();
    m_elmtV.clear();
    m_elmtA.clear();

    m_local_elmtsoln.m_gpU.clear();
    m_local_elmtsoln.m_gpUold.clear();
    m_local_elmtsoln.m_gpUolder.clear();
    m_local_elmtsoln.m_gpV.clear();
    m_local_elmtsoln.m_gpA.clear();

    m_local_elmtsoln.m_gpGradU.clear();
    m_local_elmtsoln.m_gpGradUold.clear();
    m_local_elmtsoln.m_gpGradUolder.clear();

    m_local_elmtsoln.m_gpGradV.clear();

    m_local_elmtsoln.m_gpgradu.clear();
    m_local_elmtsoln.m_gpgradv.clear();

    m_local_elmtinfo.m_gpCoords0=0.0;
    m_local_elmtinfo.m_gpCoords =0.0;

    m_nodes.clear();
    m_nodes0.clear();
}