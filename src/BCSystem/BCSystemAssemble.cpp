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
//+++ Date   : 2022.08.13
//+++ Purpose: assemble the local R and K to global one for different
//+++          boundary conditions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::assembleLocalResidual2Global(const int &Dofs,const vector<int> &DofIDs,
                                            const int &I,const double &JxW,
                                            const DofHandler &t_DofHandler,
                                            const VectorXd &LocalR,Vector &RHS){
    int iInd;
    for(int i=0;i<Dofs;i++){
        iInd=t_DofHandler.getIthNodeJthDofID(I,DofIDs[i]);
        RHS.addValue(iInd,LocalR(i+1)*JxW);
    }
}

void BCSystem::assembleLocalJacobian2Global(const int &Dofs,const vector<int> &DofIDs,
                                            const int &I,const int &J,
                                            const double &JxW,
                                            const DofHandler &t_DofHandler,
                                            const MatrixXd &LocalK,
                                            SparseMatrix &AMATRIX){
    int iInd,jInd;
    for(int i=0;i<Dofs;i++){
        iInd=t_DofHandler.getIthNodeJthDofID(I,DofIDs[i]);
        for(int j=0;j<Dofs;j++){
            jInd=t_DofHandler.getIthNodeJthDofID(J,DofIDs[j]);
            AMATRIX.addValue(iInd,jInd,LocalK(i+1,j+1)*JxW);
        }
    }
}