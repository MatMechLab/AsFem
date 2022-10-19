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
//+++ Date   : 2022.07.22
//+++ Purpose: initialize the bulk element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

void BulkElmtSystem::init(const Mesh &t_mesh){
    m_bulk_elmts=t_mesh.getBulkMeshBulkElmtsNum();
    m_elemental_elmtblock_id.clear();
    m_elemental_elmtblock_id.resize(m_bulk_elmts);
    int ee;
    for(int i=1;i<=m_elmtblock_num;i++){
        for(const auto &name:m_elmtblock_list[i-1].m_domain_namelist){
            for(int e=1;e<=t_mesh.getBulkMeshElmtsNumViaPhyName(name);e++){
                ee=t_mesh.getBulkMeshIthBulkElmtIDViaPhyName(name,e);// global element id
                m_elemental_elmtblock_id[ee-1].push_back(i);// here each bulk element could be assigned by
                                                            // multiple 'elmts' block from the input file
            }
        }
    }
    
}