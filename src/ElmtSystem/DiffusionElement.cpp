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
//+++ Date   : 2021.04.10
//+++ Purpose: implement the residual and jacobian for diffusion equation
//+++          dc/dt=div(D*grad(c))
//+++          D could be nonlinear function !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/DiffusionElement.h"

void DiffusionElement::computeAll(const FECalcType &CalcType,
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
    else{
        MessagePrinter::printErrorTxt("unsupported calculation type in DiffusionElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void DiffusionElement::computeResidual(const LocalElmtInfo &ElmtInfo,
                                       const LocalElmtSolution &ElmtSoln,
                                       const LocalShapeFun &Shp,
                                       const MaterialsContainer &MateOld,
                                       const MaterialsContainer &Mate,
                                       VectorXd &LocalR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(ElmtInfo.m_Dt||Shp.m_Test||MateOld.getScalarMaterialsNum()) {}

    LocalR(1)=ElmtSoln.m_QpV[1]*Shp.m_Test
             +Mate.ScalarMaterial("D")*(ElmtSoln.m_QpGradU[1]*Shp.m_GradTest);

}
//*****************************************************************************
void DiffusionElement::computeJacobian(const LocalElmtInfo &ElmtInfo,
                                       const double (&Ctan)[3],
                                       const LocalElmtSolution &ElmtSoln,
                                       const LocalShapeFun &Shp,
                                       const MaterialsContainer &MateOld,
                                       const MaterialsContainer &Mate,
                                       MatrixXd &LocalK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(ElmtInfo.m_Dt||MateOld.getScalarMaterialsNum()){}

    LocalK(1,1)=Shp.m_Trial*Shp.m_Test*Ctan[1]
               +Mate.ScalarMaterial("dDdc")*Shp.m_Trial*(ElmtSoln.m_QpGradU[1]*Shp.m_GradTest)*Ctan[0]
               +Mate.ScalarMaterial("D")*(Shp.m_GradTrial*Shp.m_GradTest)*Ctan[0];

}