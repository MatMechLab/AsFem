//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
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
    // TODO: we need a better assemble algorithm for complex coupling case !!!
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
//**********************************************************************
void FESystem::AssembleLocalMaterialsToGlobal(const int &e,const int &ngp,const int &gpInd,const Materials &mate,SolutionSystem &solutionSystem){
    solutionSystem._ScalarMaterials[(e-1)*ngp+gpInd-1]=mate.GetScalarMate();
    solutionSystem._VectorMaterials[(e-1)*ngp+gpInd-1]=mate.GetVectorMate();
    solutionSystem._Rank2TensorMaterials[(e-1)*ngp+gpInd-1]=mate.GetRank2Mate();
    solutionSystem._Rank4TensorMaterials[(e-1)*ngp+gpInd-1]=mate.GetRank4Mate();
}
void FESystem::AssembleLocalHistToGlobal(const int &e,const int &ngp,SolutionSystem &solutionSystem){
    for(int gpInd=0;gpInd<ngp;gpInd++){
        solutionSystem._ScalarMaterialsOld[(e-1)*ngp+gpInd]=solutionSystem._ScalarMaterials[(e-1)*ngp+gpInd];
        solutionSystem._VectorMaterialsOld[(e-1)*ngp+gpInd]=solutionSystem._VectorMaterials[(e-1)*ngp+gpInd];
        solutionSystem._Rank2TensorMaterialsOld[(e-1)*ngp+gpInd]=solutionSystem._Rank2TensorMaterials[(e-1)*ngp+gpInd];
        solutionSystem._Rank4TensorMaterialsOld[(e-1)*ngp+gpInd]=solutionSystem._Rank4TensorMaterials[(e-1)*ngp+gpInd];
    }
}
