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
//+++ Date   : 2020.11.29
//+++ Purpose: Implement the general assemble process in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/BulkFESystem.h"


//********************************************************
//*** for local to global assemble
//********************************************************
// max sub-element dofs per node (a sub-element block can never have more dofs than the
// problem's total nodal dofs, which is tiny); used to size fixed local assembly buffers.
static const int ASFEM_MAX_SUBDOFS=27;

void BulkFESystem::assembleLocalResidual2GlobalR(const int &Dofs,
                                                 const vector<int> &DofIDs,
                                                 const int &GlobalNodeID,
                                                 const DofHandler &t_DofHandler,
                                                 const double &JxW,
                                                 const VectorXd &SubR,
                                                 Vector &RHS){
    if(Dofs>ASFEM_MAX_SUBDOFS){
        MessagePrinter::printErrorTxt("sub-element dofs exceed ASFEM_MAX_SUBDOFS in residual assemble");
        MessagePrinter::exitAsFem();
    }
    // batch the per-node contributions into one VecSetValues call (indices are 0-based here,
    // unlike the per-entry addValue which subtracts 1 internally).
    int iInds[ASFEM_MAX_SUBDOFS];
    double vals[ASFEM_MAX_SUBDOFS];
    for(int i=0;i<Dofs;i++){
        iInds[i]=t_DofHandler.getIthNodeJthDofID(GlobalNodeID,DofIDs[i])-1;
        vals[i]=SubR(i+1)*JxW;
    }
    RHS.addValues(Dofs,iInds,vals);
}
void BulkFESystem::assembleLocalJacobian2GlobalK(const int &Dofs,
                                                 const vector<int> &DofIDs,
                                                 const int &GlobalNodeIDI,
                                                 const int &GlobalNodeIDJ,
                                                 const double &JxW,
                                                 const DofHandler &t_DofHandler,
                                                 const MatrixXd &SubK,
                                                 SparseMatrix &AMATRIX){
    if(Dofs>ASFEM_MAX_SUBDOFS){
        MessagePrinter::printErrorTxt("sub-element dofs exceed ASFEM_MAX_SUBDOFS in jacobian assemble");
        MessagePrinter::exitAsFem();
    }
    // batch the Dofs x Dofs block of this (node I, node J) pair into a single MatSetValues call
    // (row-major, 0-based indices) instead of one MatSetValue per entry.
    int iInds[ASFEM_MAX_SUBDOFS],jInds[ASFEM_MAX_SUBDOFS];
    double vals[ASFEM_MAX_SUBDOFS*ASFEM_MAX_SUBDOFS];
    for(int i=0;i<Dofs;i++) iInds[i]=t_DofHandler.getIthNodeJthDofID(GlobalNodeIDI,DofIDs[i])-1;
    for(int j=0;j<Dofs;j++) jInds[j]=t_DofHandler.getIthNodeJthDofID(GlobalNodeIDJ,DofIDs[j])-1;
    for(int i=0;i<Dofs;i++){
        for(int j=0;j<Dofs;j++){
            const double kij=SubK(i+1,j+1);
            vals[i*Dofs+j]=kij*JxW;
            const double ak=kij<0.0?-kij:kij;
            if(ak>m_MaxKmatCoeff) m_MaxKmatCoeff=ak;
        }
    }
    AMATRIX.addValues(Dofs,iInds,Dofs,jInds,vals);
}