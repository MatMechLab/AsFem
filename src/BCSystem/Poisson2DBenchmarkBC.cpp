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
//+++ Date   : 2022.08.26
//+++ Purpose: implement the dirichlet boundary condition for 2d
//+++          poisson benchmark test
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/Poisson2DBenchmarkBC.h"

Poisson2DBenchmarkBC::Poisson2DBenchmarkBC(){
    m_LocalU.resize(11,0.0);// the maximum dofs of each bc block is 10!
}

void Poisson2DBenchmarkBC::computeBCValue(const FECalcType &CalcType,
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
    else if (CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
        for(int i=0;i<static_cast<int>(DofIDs.size());i++){
            RHS.insertValue(DofIDs[i],0.0);
            K.insertValue(DofIDs[i],DofIDs[i],Penalty);
        }
    }

    computeU(BCValue,Params,DofIDs,ElmtInfo,ElmtSoln,m_LocalU);
    for(int i=0;i<static_cast<int>(DofIDs.size());i++){
        U.insertValue(DofIDs[i],m_LocalU(i+1));
    }

}

void Poisson2DBenchmarkBC::computeU(const double &BCValue,
                                    const nlohmann::json &Params,
                                    const vector<int> &DofIDs,
                                    const LocalElmtInfo &ElmtInfo,
                                    const LocalElmtSolution &ElmtSoln,
                                    VectorXd &LocalU){
    if(Params.size()||BCValue||ElmtSoln.m_QpU.size()){}
    double x,y;
    x=ElmtInfo.m_QpCoords0(1);
    y=ElmtInfo.m_QpCoords0(2);
    for(int i=1;i<=static_cast<int>(DofIDs.size());i++){
        LocalU(i)=1.0+x*x+2.0*y*y;
    }
}