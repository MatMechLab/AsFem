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
//+++ Date   : 2022.08.13
//+++ Purpose: assemble the local R and K to global one for different
//+++          boundary conditions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::assembleLocalResidual2Global(const int &dofs,const vector<int> &dofids,
                                            const int &I,const double &jxw,
                                            const DofHandler &dofhandler,
                                            const VectorXd &localR,Vector &RHS){
    int iInd;
    for(int i=0;i<dofs;i++){
        iInd=dofhandler.getIthNodeJthDofID(I,dofids[i]);
        RHS.addValue(iInd,localR(i+1)*jxw);
    }
}

void BCSystem::assembleLocalJacobian2Global(const int &dofs,const vector<int> &dofids,
                                            const int &I,const int &J,
                                            const double &jxw,
                                            const DofHandler &dofhandler,
                                            const MatrixXd &localK,
                                            SparseMatrix &AMATRIX){
    int iInd,jInd;
    for(int i=0;i<dofs;i++){
        iInd=dofhandler.getIthNodeJthDofID(I,dofids[i]);
        for(int j=0;j<dofs;j++){
            jInd=dofhandler.getIthNodeJthDofID(J,dofids[j]);
            AMATRIX.addValue(iInd,jInd,localK(i+1,j+1)*jxw);
        }
    }
}