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
//+++ Purpose: Implement the general assemble process in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/BulkFESystem.h"


//********************************************************
//*** for local to global assemble
//********************************************************
void BulkFESystem::assembleLocalResidual2GlobalR(const int &t_dofs,const vector<int> &t_dofsid,
                                       const int &t_globalnodeid,
                                       const DofHandler &t_dofhandler,
                                       const double &jxw,
                                       const VectorXd &t_subR,
                                       Vector &RHS){
    int iInd;
    for(int i=0;i<t_dofs;i++){
        iInd=t_dofhandler.getIthNodeJthDofID(t_globalnodeid,t_dofsid[i]);
        RHS.addValue(iInd,t_subR(i+1)*jxw);
    }
}
void BulkFESystem::assembleLocalJacobian2GlobalK(const int &t_dofs,const vector<int> &t_dofsid,
                                       const int &t_globalnodeidI,const int &t_globalnodeidJ,
                                       const double &jxw,
                                       const DofHandler &t_dofhandler,
                                       const MatrixXd &t_subK,
                                       SparseMatrix &AMATRIX){
    int iInd,jInd;
    for(int i=0;i<t_dofs;i++){
        iInd=t_dofhandler.getIthNodeJthDofID(t_globalnodeidI,t_dofsid[i]);
        for(int j=0;j<t_dofs;j++){
            jInd=t_dofhandler.getIthNodeJthDofID(t_globalnodeidJ,t_dofsid[j]);
            AMATRIX.addValue(iInd,jInd,t_subK(i+1,j+1)*jxw*1.0);
            if(abs(t_subK(i+1,j+1))>m_max_k_coeff) m_max_k_coeff=abs(t_subK(i+1,j+1));
        }
    }
}