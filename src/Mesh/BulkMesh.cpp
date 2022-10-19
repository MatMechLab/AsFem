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
//+++ Date   : 2022.05.06
//+++ Purpose: the bulk mesh class of AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/BulkMesh.h"

BulkMesh::BulkMesh(){
}
BulkMesh::BulkMesh(const BulkMesh &mesh){
    m_meshdata=mesh.m_meshdata;
}
void BulkMesh::releaseMemory(){
    m_meshdata.m_nodecoords0.clear();
    m_meshdata.m_nodecoords.clear();
    m_meshdata.m_bulkelmt_connectivity.clear();
    m_meshdata.m_bulkelmt_volume.clear();
    m_meshdata.m_pointelmt_connectivity.clear();
    m_meshdata.m_pointelmt_volume.clear();
    m_meshdata.m_lineelmt_connectivity.clear();
    m_meshdata.m_lineelmt_volume.clear();
    m_meshdata.m_surfaceelmt_connectivity.clear();
    m_meshdata.m_surfaceelmt_volume.clear();

    m_meshdata.m_phygroup_dimvec.clear();
    m_meshdata.m_phygroup_phynamevec.clear();
    m_meshdata.m_phygroup_phyidvec.clear();
    m_meshdata.m_phygroup_elmtnumvec.clear();
    m_meshdata.m_phygroup_nodesnumperelmtvec.clear();
    m_meshdata.m_phygroup_name2elmtconnvec.clear();
    m_meshdata.m_phygroup_name2bulkelmtidvec.clear();
    m_meshdata.m_phygroup_name2dimvec.clear();
    m_meshdata.m_phygroup_name2phyidvec.clear();
    m_meshdata.m_phygroup_phyid2namevec.clear();

    m_meshdata.m_nodephygroup_name2nodeidvec.clear();
    m_meshdata.m_nodephygroup_name2phyidvec.clear();
    m_meshdata.m_nodephygroup_phyid2namevec.clear();
    m_meshdata.m_nodephygroup_phynamevec.clear();
    m_meshdata.m_nodephygroup_phyidvec.clear();

}