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
//+++ Date   : 2021.10.06
//+++ Purpose: implement the dirichlet boundary condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/DirichletBC.h"

DirichletBC::DirichletBC(){
    m_LocalU.resize(11,0.0);// the maximum dofs of each bc block is 10!
}

void DirichletBC::computeBCValue(const FECalcType &CalcType,
                                 const double &Penalty,
                                 const double &BCValue,
                                 const nlohmann::json &Params,
                                 const LocalElmtInfo &ElmtInfo,
                                 const LocalElmtSolution &ElmtSoln,
                                 const vector<int> &DofIDs,
                                 Vector &U,
                                 SparseMatrix &K,
                                 Vector &RHS){
    if(CalcType==FECalcType::COMPUTERESIDUAL){
        for(int i=0;i<static_cast<int>(DofIDs.size());i++){
            RHS.insertValue(DofIDs[i],0.0);
        }
    }
    else if(CalcType==FECalcType::COMPUTEJACOBIAN){
        for(int i=0;i<static_cast<int>(DofIDs.size());i++){
            K.insertValue(DofIDs[i],DofIDs[i],Penalty);
        }
    }

    computeU(BCValue,Params,DofIDs,ElmtInfo,ElmtSoln,m_LocalU);
    for(int i=0;i<static_cast<int>(DofIDs.size());i++){
        U.insertValue(DofIDs[i],m_LocalU(i+1));
    }

}

void DirichletBC::computeU(const double &BCValue,
                           const nlohmann::json &Params,
                           const vector<int> &DofIDs,
                           const LocalElmtInfo &ElmtInfo,
                           const LocalElmtSolution &ElmtSoln,
                           VectorXd &LocalU){
    if(Params.size()||ElmtInfo.m_T||ElmtSoln.m_QpU.size()){}
    for(int i=1;i<=static_cast<int>(DofIDs.size());i++){
        LocalU(i)=BCValue;
    }
}