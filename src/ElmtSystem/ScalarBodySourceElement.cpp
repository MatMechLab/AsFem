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
//+++ Date   : 2022.11.20
//+++ Purpose: implement the residual and jacobian for scalar body
//+++          source contribution.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/ScalarBodySourceElement.h"

void ScalarBodySourceElement::computeAll(const FECalcType &CalcType,
                                const LocalElmtInfo &ElmtInfo,
                                const double (&Ctan)[3],
                                const LocalElmtSolution &ElmtSoln,
                                const LocalShapeFun &Shp,
                                const MaterialsContainer &MateOld,
                                const MaterialsContainer &Mate,
                                MatrixXd &LocalK,VectorXd &LocalR) {
    if(CalcType==FECalcType::COMPUTERESIDUAL){
        computeResidual(ElmtInfo,ElmtSoln,Shp,MateOld,Mate,LocalR);
    }
    else if(CalcType==FECalcType::COMPUTEJACOBIAN){
        computeJacobian(ElmtInfo,Ctan,ElmtSoln,Shp,MateOld,Mate,LocalK);
    }
    else if(CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN){
        computeResidual(ElmtInfo,ElmtSoln,Shp,MateOld,Mate,LocalR);
        computeJacobian(ElmtInfo,Ctan,ElmtSoln,Shp,MateOld,Mate,LocalK);
    }
    else{
        MessagePrinter::printErrorTxt("unsupported calculation type in ScalarBodySourceElement, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void ScalarBodySourceElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_Dt||soln.m_QpU[0]||shp.m_Test||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()) {}

    localR(1)=mate.ScalarMaterial("f")*shp.m_Test;

}
//*****************************************************************************
void ScalarBodySourceElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&Ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_Dt||Ctan[0]||soln.m_QpU[0]||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()){}

    localK(1,1)=mate.ScalarMaterial("dfdu")*shp.m_Trial*shp.m_Test*Ctan[0];

}