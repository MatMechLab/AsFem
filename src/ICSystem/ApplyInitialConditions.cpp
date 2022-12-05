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
//+++ Date   : 2020.07.10
//+++ Purpose: Apply the different initial conditions according to
//+++          the input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

void ICSystem::applyInitialConditions(const Mesh &t_mesh,const DofHandler &t_dofhandler,Vector &U0){
    int e,i,j,k,iInd,dofs,dim;
    int rankne,eStart,eEnd,nElmts,nNodesPerElmt;
    double icvalue;
    Vector3d nodecoords0;

    MPI_Comm_size(PETSC_COMM_WORLD,&m_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);

    PetscRandomCreate(PETSC_COMM_WORLD,&m_rnd);
    PetscRandomSetType(m_rnd,PETSCRAND);

    m_minval=0.0;m_maxval=1.0;

    for(const auto &it:m_icblock_list){
        icvalue=it.m_icvalue;
        dofs=static_cast<int>(it.m_dofIDs.size());
        if(it.m_icType==ICType::RANDOMIC){
            m_minval=JsonUtils::getValue(it.m_json_params,"minval");
            m_maxval=JsonUtils::getValue(it.m_json_params,"maxval");
            PetscRandomSetInterval(m_rnd,m_minval,m_maxval);
        }
        for(const auto &name:it.m_domainNameList){
            nElmts=t_mesh.getBulkMeshElmtsNumViaPhyName(name);
            rankne=nElmts/m_size;
            eStart=m_rank*rankne;
            eEnd=(m_rank+1)*rankne;
            if(m_rank==m_size-1) eEnd=nElmts;
            dim=t_mesh.getBulkMeshElmtDimViaPhyName(name);
            for(e=eStart;e<eEnd;e++){
                nNodesPerElmt=t_mesh.getBulkMeshIthElmtNodesNumViaPhyName(name,e+1);
                for(i=1;i<=nNodesPerElmt;i++){
                    if(it.m_icType==ICType::RANDOMIC) PetscRandomGetValue(m_rnd,&icvalue);
                    j=t_mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);// global id
                    nodecoords0(1)=t_mesh.getBulkMeshIthNodeJthCoord0(j,1);
                    nodecoords0(2)=t_mesh.getBulkMeshIthNodeJthCoord0(j,2);
                    nodecoords0(3)=t_mesh.getBulkMeshIthNodeJthCoord0(j,3);
                    runICLibs(it.m_icType,it.m_json_params,icvalue,dim,dofs,nodecoords0,m_localU);
                    for(k=1;k<=dofs;k++){
                        iInd=t_dofhandler.getIthNodeJthDofID(j,it.m_dofIDs[k-1]);
                        U0.insertValue(iInd,m_localU(k));
                    }
                }
            }
        }
    }

    U0.assemble();

    PetscRandomDestroy(&m_rnd);

}