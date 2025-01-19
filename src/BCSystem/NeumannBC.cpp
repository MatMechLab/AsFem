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
//+++ Purpose: implement the Neumann type boundary condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/NeumannBC.h"

void NeumannBC::computeBCValue(const FECalcType &CalcType,
                               const double &BCValue,
                               const nlohmann::json Params,
                               const LocalElmtInfo &ElmtInfo,
                               const LocalElmtSolution &ElmtSoln,
                               const Vector3d &Normal,
                               const LocalShapeFun &Shp,
                               const double (&Ctan)[3],
                               MatrixXd &LocalK,
                               VectorXd &LocalR){
    if(CalcType==FECalcType::COMPUTERESIDUAL){
        computeResidual(BCValue,Params,ElmtInfo,ElmtSoln,Normal,Shp,LocalR);
    }
    else if(CalcType==FECalcType::COMPUTEJACOBIAN){
        computeJacobian(BCValue,Params,ElmtInfo,ElmtSoln,Normal,Shp,Ctan,LocalK);
    }
    else if (CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
        computeResidual(BCValue,Params,ElmtInfo,ElmtSoln,Normal,Shp,LocalR);
        computeJacobian(BCValue,Params,ElmtInfo,ElmtSoln,Normal,Shp,Ctan,LocalK);
    }
    else{
        MessagePrinter::printErrorTxt("Unsupported calculation type in NeumannBC class, please check your code");
        MessagePrinter::exitAsFem();
    }
}


void NeumannBC::computeResidual(const double &BCValue,
                                const nlohmann::json Params,
                                const LocalElmtInfo &ElmtInfo,
                                const LocalElmtSolution &ElmtSoln,
                                const Vector3d &Normal,
                                const LocalShapeFun &Shp,
                                VectorXd &LocalR){
    // get rid of unused warnings
    if(Params.size()||ElmtInfo.m_Dt||ElmtSoln.m_QpU[0]||Normal(1)){}

    LocalR(1)=BCValue*Shp.m_Test;
}

void NeumannBC::computeJacobian(const double &BCValue,
                                const nlohmann::json &Params,
                                const LocalElmtInfo &ElmtInfo,
                                const LocalElmtSolution &ElmtSoln,
                                const Vector3d &Normal,
                                const LocalShapeFun &Shp,
                                const double (&Ctan)[3],
                                MatrixXd &LocalK){
    // get rid of unused warnings
    if(BCValue||Params.size()||ElmtInfo.m_Dt||ElmtSoln.m_QpU[0]||Normal(1)||Shp.m_Test){}
    LocalK(1,1)=0.0*Ctan[0]; // since bcvalue is a constant, its derivative should be zero!

}