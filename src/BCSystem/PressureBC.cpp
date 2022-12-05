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
//+++ Date   : 2022.08.25
//+++ Purpose: implement the pressure type boundary condition for
//+++          solid mechanics problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/PressureBC.h"

void PressureBC::computeBCValue(const FECalcType &calctype,const double &bcvalue,
                               const nlohmann::json parameters,
                               const LocalElmtInfo &elmtinfo,
                               const LocalElmtSolution &elmtsoln,
                               const Vector3d &normal,
                               const LocalShapeFun &shp,
                               const double (&ctan)[3],
                               MatrixXd &localK,
                               VectorXd &localR){
    if(calctype==FECalcType::COMPUTERESIDUAL){
        computeResidual(bcvalue,parameters,elmtinfo,elmtsoln,normal,shp,localR);
    }
    else if(calctype==FECalcType::COMPUTEJACOBIAN){
        computeJacobian(bcvalue,parameters,elmtinfo,elmtsoln,normal,shp,ctan,localK);
    }
    else{
        MessagePrinter::printErrorTxt("Unsupported calculation type in PressureBC class, please check your code");
        MessagePrinter::exitAsFem();
    }
}


void PressureBC::computeResidual(const double &bcvalue,
                                const nlohmann::json parameters,
                                const LocalElmtInfo &elmtinfo,
                                const LocalElmtSolution &elmtsoln,
                                const Vector3d &normal,
                                const LocalShapeFun &shp,
                                VectorXd &localR){
    // get rid of unused warnings
    if(bcvalue||elmtinfo.m_dim||elmtsoln.m_gpU[0]){}

    m_pressure=JsonUtils::getValue(parameters,"pressure");
    m_component=static_cast<int>(JsonUtils::getValue(parameters,"component"));
    if(m_component<1||m_component>3){
        MessagePrinter::printErrorTxt("Invalid component(="+to_string(m_component)+") for PressureBC, please check your input file");
        MessagePrinter::exitAsFem();
    }
    m_traction=normal*m_pressure;

    localR(1)=m_traction(m_component)*shp.m_test;
}

void PressureBC::computeJacobian(const double &bcvalue,
                                const nlohmann::json &parameters,
                                const LocalElmtInfo &elmtinfo,
                                const LocalElmtSolution &elmtsoln,
                                const Vector3d &normal,
                                const LocalShapeFun &shp,
                                const double (&ctan)[3],
                                MatrixXd &localK){
    // get rid of unused warnings
    if(bcvalue||parameters.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||normal(1)||shp.m_test||ctan[0]){}
    localK(1,1)=0.0; // since bcvalue is a constant, its derivative should be zero!

}