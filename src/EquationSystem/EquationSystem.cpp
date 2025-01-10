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
//+++ Date   : 2022.07.24
//+++ Purpose: the equation system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "EquationSystem/EquationSystem.h"

EquationSystem::EquationSystem(){
    m_Dofs=0;
    m_Allocated=false;
}

void EquationSystem::init(const DofHandler &t_dofHandler){
    m_Dofs=t_dofHandler.getActiveDofs();
    m_RHS.resize(m_Dofs,0.0);
    m_AMATRIX.resize(m_Dofs,m_Dofs,t_dofHandler.getMaxRowNNZ());
    m_Allocated=true;
}

void EquationSystem::createSparsityPattern(const DofHandler &t_dofHandler){
    vector<int> eldofs;
    vector<double> elvals;
    int n;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    eldofs.resize(t_dofHandler.getMaxDofsPerElmt()+1,0);
    elvals.resize(t_dofHandler.getMaxDofsPerElmt()+1,0.0);

    n=t_dofHandler.getMaxDofsPerElmt();
    for(int e=1;e<=t_dofHandler.getLocalBulkElmtsNum();e++){
        t_dofHandler.getIthLocalBulkElmtDofIDs(e,eldofs);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                m_AMATRIX.addValue(eldofs[i],eldofs[j],elvals[i]);
            }
        }
    }
    m_AMATRIX.assemble();
    m_AMATRIX.disableReallocation();// the following operation can not modify the sparsity pattern anymore!!!

    eldofs.clear();
    elvals.clear();
}

void EquationSystem::releaseMemory(){
    m_RHS.releaseMemory();
    m_AMATRIX.releaseMemory();
}