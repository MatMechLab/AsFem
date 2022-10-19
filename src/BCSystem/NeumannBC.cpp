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
//+++ Date   : 2021.10.06
//+++ Purpose: implement the Neumann type boundary condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/NeumannBC.h"

void NeumannBC::computeBCValue(const FECalcType &calctype,const double &bcvalue,
                               const nlohmann::json json,
                               const LocalElmtInfo &elmtinfo,
                               const LocalElmtSolution &elmtsoln,
                               const Vector3d &normal,
                               const LocalShapeFun &shp,
                               const double (&ctan)[3],
                               MatrixXd &localK,
                               VectorXd &localR){
    if(calctype==FECalcType::COMPUTERESIDUAL){
        computeResidual(bcvalue,json,elmtinfo,elmtsoln,normal,shp,localR);
    }
    else if(calctype==FECalcType::COMPUTEJACOBIAN){
        computeJacobian(bcvalue,json,elmtinfo,elmtsoln,normal,shp,ctan,localK);
    }
    else{
        MessagePrinter::printErrorTxt("Unsupported calculation type in NeumannBC class, please check your code");
        MessagePrinter::exitAsFem();
    }
}


void NeumannBC::computeResidual(const double &bcvalue,
                                const nlohmann::json json,
                                const LocalElmtInfo &elmtinfo,
                                const LocalElmtSolution &elmtsoln,
                                const Vector3d &normal,
                                const LocalShapeFun &shp,
                                VectorXd &localR){
    // get rid of unused warnings
    if(json.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||normal(1)){}

    localR(1)=bcvalue*shp.m_test;
}

void NeumannBC::computeJacobian(const double &bcvalue,
                                const nlohmann::json &json,
                                const LocalElmtInfo &elmtinfo,
                                const LocalElmtSolution &elmtsoln,
                                const Vector3d &normal,
                                const LocalShapeFun &shp,
                                const double (&ctan)[3],
                                MatrixXd &localK){
    // get rid of unused warnings
    if(bcvalue||json.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||normal(1)||shp.m_test){}
    localK(1,1)=0.0*ctan[0]; // since bcvalue is a constant, its derivative should be zero!

}