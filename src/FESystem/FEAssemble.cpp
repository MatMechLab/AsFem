//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.11.30
//+++ Purpose: Implement the different assemble functions for different
//+++          cases, i.e., residual assemble and jacobian assemble
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/FESystem.h"

void FESystem::AssembleSubResidualToLocalResidual(const int &ndofspernode,const int &dofs,const int &iInd,
                                            const VectorXd &subR,VectorXd &localR){
    for(int i=1;i<=dofs;i++){
        localR((iInd-1)*ndofspernode+i)+=subR(i);
    }
}
//***********************************
void FESystem::AccumulateLocalResidual(const int &dofs,const vector<double> &dofsactiveflag,const double &JxW,
                                const VectorXd &localR,vector<double> &sumR){
    for(int i=1;i<=dofs;i++){
        sumR[i-1]+=localR(i)*JxW*dofsactiveflag[i-1];
    }
}
//*****************************************
void FESystem::AssembleLocalResidualToGlobalResidual(const int &ndofs,const vector<int> &dofindex,
                                            const vector<double> &residual,Vec &rhs){
    VecSetValues(rhs,ndofs,dofindex.data(),residual.data(),ADD_VALUES);
}
//*************************************************************
void FESystem::AssembleSubJacobianToLocalJacobian(const int &ndofspernode,
                                            const int &iInd,const int &jInd,
                                            const MatrixXd &subK,MatrixXd &localK){
    for(int i=1;i<=ndofspernode;i++){
        for(int j=1;j<=ndofspernode;j++){
            localK((iInd-1)*ndofspernode+i,(jInd-1)*ndofspernode+j)+=subK(i,j);
        }
    }
}
void FESystem::AccumulateLocalJacobian(const int &dofs,const vector<double> &dofsactiveflag,const double &JxW,
                                const MatrixXd &localK,vector<double> &sumK){
    for(int i=1;i<=dofs;i++){
        if(dofsactiveflag[i-1]>0.0){
            for(int j=1;j<=dofs;j++){
                sumK[(i-1)*dofs+j-1]+=localK(i,j)*JxW;
                if(localK(i,j)*JxW>_MaxKMatrixValue) _MaxKMatrixValue=localK(i,j)*JxW;
            }
        }
    }
}
void FESystem::AssembleLocalJacobianToGlobalJacobian(const int &ndofs,const vector<int> &dofindex,
                                            const vector<double> &jacobian,Mat &K){
    MatSetValues(K,ndofs,dofindex.data(),ndofs,dofindex.data(),jacobian.data(),ADD_VALUES);
}