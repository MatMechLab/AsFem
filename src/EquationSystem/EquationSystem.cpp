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
//+++ Date   : 2022.07.24
//+++ Purpose: the equation system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "EquationSystem/EquationSystem.h"

EquationSystem::EquationSystem(){
    m_dofs=0;
    m_allocated=false;
}

void EquationSystem::init(const DofHandler &t_dofHandler){
    m_dofs=t_dofHandler.getActiveDofs();
    m_rhs.resize(m_dofs,0.0);
    m_amatrix.resize(m_dofs,m_dofs,t_dofHandler.getMaxNNZ());
    m_allocated=true;
}

void EquationSystem::createSparsityPattern(const DofHandler &t_dofHandler){
    vector<int> eldofs;
    vector<double> elvals;
    int n;

    PetscMPIInt rank,size;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    eldofs.resize(t_dofHandler.getMaxDofsPerElmt()+1,0);
    elvals.resize(t_dofHandler.getMaxDofsPerElmt()+1,0.0);

    int rankne=t_dofHandler.getBulkElmtsNum()/size;
    int eStart=rank*rankne;
    int eEnd=(rank+1)*rankne;
    if(rank==size-1) eEnd=t_dofHandler.getBulkElmtsNum();
    int e;
    for(int ee=eStart;ee<eEnd;++ee){
        e=ee+1;
        t_dofHandler.getIthBulkElmtDofIDs0(e,eldofs);// index start from 0
        n=t_dofHandler.getIthBulkElmtDofsNum(e);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                m_amatrix.addValue(eldofs[i]+1,eldofs[j]+1,elvals[i]);
            }
        }
    }
    m_amatrix.assemble();
    m_amatrix.disableReallocation();// the following operation can not modify the sparsity pattern anymore!!!

    eldofs.clear();
    elvals.clear();
}

void EquationSystem::releaseMemory(){
    m_rhs.releaseMemory();
    m_amatrix.releaseMemory();
}