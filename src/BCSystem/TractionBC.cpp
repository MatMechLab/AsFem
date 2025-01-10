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
//+++ Date   : 2022.08.25
//+++ Purpose: implement the traction boundary condition for
//+++          solid mechanics problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/TractionBC.h"

void TractionBC::computeBCValue(const FECalcType &CalcType,const double &BCValue,
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
    else{
        MessagePrinter::printErrorTxt("Unsupported calculation type in TractionBC class, please check your code");
        MessagePrinter::exitAsFem();
    }
}


void TractionBC::computeResidual(const double &BCValue,
                                 const nlohmann::json Params,
                                 const LocalElmtInfo &ElmtInfo,
                                 const LocalElmtSolution &ElmtSoln,
                                 const Vector3d &Normal,
                                 const LocalShapeFun &Shp,
                                 VectorXd &LocalR){
    // get rid of unused warnings
    if(BCValue||ElmtInfo.m_Dim||Normal(1)||ElmtSoln.m_QpU[0]){}

    m_traction=JsonUtils::getVector(Params,"traction");
    m_component=static_cast<int>(JsonUtils::getValue(Params,"component"));
    if(m_component<1||m_component>3){
        MessagePrinter::printErrorTxt("Invalid component(="+to_string(m_component)+") for TractionBC, please check your input file");
        MessagePrinter::exitAsFem();
    }
    LocalR(1)=m_traction(m_component)*-1.0*Shp.m_Test;
}

void TractionBC::computeJacobian(const double &BCValue,
                                 const nlohmann::json &Params,
                                 const LocalElmtInfo &ElmtInfo,
                                 const LocalElmtSolution &ElmtSoln,
                                 const Vector3d &Normal,
                                 const LocalShapeFun &Shp,
                                 const double (&Ctan)[3],
                                 MatrixXd &LocalK){
    // get rid of unused warnings
    if(BCValue||Params.size()||ElmtInfo.m_Dt||ElmtSoln.m_QpU[0]||Normal(1)||Shp.m_Test||Ctan[0]){}
    LocalK(1,1)=0.0; // since bcvalue is a constant, its derivative should be zero!

}