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
//+++ Purpose: implement the residual and jacobian for cahn-hilliard 
//+++          equation
//+++          dc/dt=div(M*grad(u))
//+++             mu=deltaF/deltaC
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/CahnHilliardElement.h"

void CahnHilliardElement::computeAll(const FECalcType &CalcType,
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
        MessagePrinter::printErrorTxt("unsupported calculation type in CahnHilliardElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void CahnHilliardElement::computeResidual(const LocalElmtInfo &ElmtInfo,
                                          const LocalElmtSolution &ElmtSoln,
                                          const LocalShapeFun &Shp,
                                          const MaterialsContainer &MateOld,
                                          const MaterialsContainer &Mate,
                                          VectorXd &LocalR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(ElmtInfo.m_Dt||MateOld.getScalarMaterialsNum()) {}
    // R_c
    LocalR(1)=ElmtSoln.m_QpV[1]*Shp.m_Test
             +Mate.ScalarMaterial("M")*(ElmtSoln.m_QpGradU[2]*Shp.m_GradTest);
    // R_mu
    LocalR(2)=ElmtSoln.m_QpU[2]*Shp.m_Test
             -Mate.ScalarMaterial("dFdC")*Shp.m_Test
             -Mate.ScalarMaterial("kappa")*ElmtSoln.m_QpGradU[1]*Shp.m_GradTest;

}
//*****************************************************************************
void CahnHilliardElement::computeJacobian(const LocalElmtInfo &ElmtInfo,
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
    // K_c,c
    LocalK(1,1)=Shp.m_Trial*Shp.m_Test*Ctan[1]
               +Mate.ScalarMaterial("dMdC")*Shp.m_Trial*(ElmtSoln.m_QpGradU[2]*Shp.m_GradTest)*Ctan[0];
    // K_c,mu
    LocalK(1,2)=Mate.ScalarMaterial("M")*(Shp.m_GradTrial*Shp.m_GradTest)*Ctan[0];

    // K_mu,c
    LocalK(2,1)=-Mate.ScalarMaterial("d2FdC2")*Shp.m_Trial*Shp.m_Test*Ctan[0]
                -Mate.ScalarMaterial("kappa")*Shp.m_GradTrial*Shp.m_GradTest*Ctan[0];
    // K_mu,mu
    LocalK(2,2)=Shp.m_Trial*Shp.m_Test*Ctan[0];

}