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
//+++ Date   : 2020.07.10
//+++ Purpose: Apply the different initial conditions according to
//+++          the input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

#include <cstdint>

namespace{
    /**
     * Deterministic, partition-independent pseudo-random value in [0,1) derived purely
     * from a global id (splitmix64 mix). Seeding the random IC by the GLOBAL node id makes
     * the random field identical regardless of how the mesh is split across mpi ranks, and
     * idempotent for shared/ghost nodes (the same node always gets the same value). This
     * replaces the previous PetscRandom draw, which depended on the local cell-node visiting
     * order and therefore produced different fields for different cpu counts.
     * @param id the global node id used as the seed
     */
    inline double randomICUniformFromGlobalID(uint64_t id){
        uint64_t x=id+0x9E3779B97F4A7C15ULL;
        x=(x^(x>>30))*0xBF58476D1CE4E5B9ULL;
        x=(x^(x>>27))*0x94D049BB133111EBULL;
        x= x^(x>>31);
        return static_cast<double>(x>>11)*(1.0/9007199254740992.0);// 53-bit mantissa / 2^53
    }
}

void ICSystem::applyInitialConditions(const FECell &t_fecell,const DofHandler &t_dofhandler,Vector &U0){
    int iInd,GlobalNodeID,dofs,dim;
    int nNodesPerElmt;
    double icvalue;
    Vector3d nodecoords0;
    iInd=dofs=dim=1;
    nNodesPerElmt=1;
    icvalue=1.0;

    m_minval=0.0;m_maxval=1.0;

    for(const auto &block:m_icblock_list){
        icvalue=block.m_ICValue;
        dofs=static_cast<int>(block.m_DofIDs.size());
        if(block.m_ICType==ICType::RANDOMIC){
            m_minval=JsonUtils::getValue(block.m_Params,"minval");
            m_maxval=JsonUtils::getValue(block.m_Params,"maxval");
        }
        for(const auto &name:block.m_DomainNameList){
            for (const auto &cell:t_fecell.getLocalMeshCellVectorCopyViaPhyName(name)) {
                nNodesPerElmt=cell.NodesNumPerElmt;
                dim=cell.Dim;
                for (int i=1;i<=nNodesPerElmt;i++) {
                    GlobalNodeID=cell.ElmtConn[i-1];
                    if (block.m_ICType==ICType::RANDOMIC) {
                        // seed by the GLOBAL node id so the random field is reproducible
                        // across cpu counts and identical for shared/ghost nodes
                        icvalue=m_minval+(m_maxval-m_minval)*randomICUniformFromGlobalID(static_cast<uint64_t>(GlobalNodeID));
                    }
                    nodecoords0(1)=cell.ElmtNodeCoords0(i,1);
                    nodecoords0(2)=cell.ElmtNodeCoords0(i,2);
                    nodecoords0(3)=cell.ElmtNodeCoords0(i,3);
                    runICLibs(block.m_ICType,block.m_Params,icvalue,dim,dofs,nodecoords0,m_localU);
                    for (int k=1;k<=dofs;k++) {
                        iInd=t_dofhandler.getIthNodeJthDofID(GlobalNodeID,block.m_DofIDs[k-1]);
                        U0.insertValue(iInd,m_localU(k));
                    }
                }
            }
        }
    }

    U0.assemble();

}