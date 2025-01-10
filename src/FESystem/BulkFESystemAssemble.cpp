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
void BulkFESystem::assembleLocalResidual2GlobalR(const int &Dofs,
                                                 const vector<int> &DofIDs,
                                                 const int &GlobalNodeID,
                                                 const DofHandler &t_DofHandler,
                                                 const double &JxW,
                                                 const VectorXd &SubR,
                                                 Vector &RHS){
    int iInd;
    for(int i=0;i<Dofs;i++){
        iInd=t_DofHandler.getIthNodeJthDofID(GlobalNodeID,DofIDs[i]);
        RHS.addValue(iInd,SubR(i+1)*JxW);
    }
}
void BulkFESystem::assembleLocalJacobian2GlobalK(const int &Dofs,
                                                 const vector<int> &DofIDs,
                                                 const int &GlobalNodeIDI,
                                                 const int &GlobalNodeIDJ,
                                                 const double &JxW,
                                                 const DofHandler &t_DofHandler,
                                                 const MatrixXd &SubK,
                                                 SparseMatrix &AMATRIX){
    int iInd,jInd;
    for(int i=0;i<Dofs;i++){
        iInd=t_DofHandler.getIthNodeJthDofID(GlobalNodeIDI,DofIDs[i]);
        for(int j=0;j<Dofs;j++){
            jInd=t_DofHandler.getIthNodeJthDofID(GlobalNodeIDJ,DofIDs[j]);
            AMATRIX.addValue(iInd,jInd,SubK(i+1,j+1)*JxW*1.0);
            if(abs(SubK(i+1,j+1))>m_MaxKmatCoeff) m_MaxKmatCoeff=abs(SubK(i+1,j+1));
        }
    }
}