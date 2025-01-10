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

void ICSystem::applyInitialConditions(const FECell &t_fecell,const DofHandler &t_dofhandler,Vector &U0){
    int iInd,GlobalNodeID,dofs,dim;
    int nNodesPerElmt;
    double icvalue;
    Vector3d nodecoords0;
    iInd=dofs=dim=1;
    nNodesPerElmt=1;
    icvalue=1.0;

    PetscRandomCreate(PETSC_COMM_WORLD,&m_rnd);
    PetscRandomSetType(m_rnd,PETSCRAND);

    m_minval=0.0;m_maxval=1.0;

    for(const auto &block:m_icblock_list){
        icvalue=block.m_ICValue;
        dofs=static_cast<int>(block.m_DofIDs.size());
        if(block.m_ICType==ICType::RANDOMIC){
            m_minval=JsonUtils::getValue(block.m_Params,"minval");
            m_maxval=JsonUtils::getValue(block.m_Params,"maxval");
            PetscRandomSetInterval(m_rnd,m_minval,m_maxval);
        }
        for(const auto &name:block.m_DomainNameList){
            for (const auto &cell:t_fecell.getLocalMeshCellVectorCopyViaPhyName(name)) {
                nNodesPerElmt=cell.NodesNumPerElmt;
                dim=cell.Dim;
                for (int i=1;i<=nNodesPerElmt;i++) {
                    if (block.m_ICType==ICType::RANDOMIC) {
                        PetscRandomGetValue(m_rnd,&icvalue);
                    }
                    GlobalNodeID=cell.ElmtConn[i-1];
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

    PetscRandomDestroy(&m_rnd);

}