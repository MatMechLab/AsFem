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
//+++ Date   : 2022.11.12
//+++ Purpose: Calculate the residual and jacobian for dendrite growth model
//+++          based on Kobayashi's paper
//+++ Ref    : Modeling and numerical simulations of dendritic crystal growth
//+++ DOI    : https://doi.org/10.1016/0167-2789(93)90120-P
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/KobayashiElement.h"

void KobayashiElement::computeAll(const FECalcType &CalcType,
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
        MessagePrinter::printErrorTxt("unsupported calculation type in KobayashiElement, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void KobayashiElement::computeResidual(const LocalElmtInfo &ElmtInfo,
                                       const LocalElmtSolution &ElmtSoln,
                                       const LocalShapeFun &Shp,
                                       const MaterialsContainer &MateOld,
                                       const MaterialsContainer &Mate,
                                       VectorXd &LocalR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(MateOld.getScalarMaterialsNum()) {}

    if(ElmtInfo.m_Dim!=2){
        MessagePrinter::printErrorTxt("Kobayashi element only works for 2d case");
        MessagePrinter::exitAsFem();
    }

    L=Mate.ScalarMaterial("L");
    K=Mate.ScalarMaterial("K");
    dK=Mate.ScalarMaterial("dK");
    dFdeta=Mate.ScalarMaterial("dFdeta");
    Latent=Mate.ScalarMaterial("Latent-heat");
    // For R_eta
    V(1)=-ElmtSoln.m_QpGradU[1](2);
    V(2)= ElmtSoln.m_QpGradU[1](1);
    V(3)= 0.0;
    LocalR(1)=ElmtSoln.m_QpV[1]*Shp.m_Test
             +L*K*dK*(V*Shp.m_GradTest)
             +L*K*K*(ElmtSoln.m_QpGradU[1]*Shp.m_GradTest)
             +L*dFdeta*Shp.m_Test;
    // For R_T
    LocalR(2)=ElmtSoln.m_QpV[2]*Shp.m_Test
             +ElmtSoln.m_QpGradU[2]*Shp.m_GradTest
             -Latent*ElmtSoln.m_QpV[1]*Shp.m_Test;

}
//*****************************************************************************
void KobayashiElement::computeJacobian(const LocalElmtInfo &ElmtInfo,
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
    L=Mate.ScalarMaterial("L");
    K=Mate.ScalarMaterial("K");
    dK=Mate.ScalarMaterial("dK");
    Latent=Mate.ScalarMaterial("Latent-heat");
    d2Fdeta2=Mate.ScalarMaterial("d2Fdeta2");
    d2FdetadT=Mate.ScalarMaterial("d2FdetadT");
    dKdGradEta=Mate.VectorMaterial("dKdGradEta");
    ddKdGradEta=Mate.VectorMaterial("ddKdGradEta");
    //************************************************
    V(1)=-ElmtSoln.m_QpGradU[1](2);
    V(2)= ElmtSoln.m_QpGradU[1](1);
    V(3)= 0.0;
    //************************
    dV(1)=-Shp.m_GradTrial(2);
    dV(2)= Shp.m_GradTrial(1);
    dV(3)= 0.0;
    // K_eta,eta
    LocalK(1,1)=Shp.m_Trial*Shp.m_Test*Ctan[1]
               +L*(dKdGradEta*Shp.m_GradTrial)*dK*(V*Shp.m_GradTest)*Ctan[0]
               +L*K*(ddKdGradEta*Shp.m_GradTrial)*(V*Shp.m_GradTest)*Ctan[0]
               +L*K*dK*(dV*Shp.m_GradTest)*Ctan[0]
               +L*2*K*(dKdGradEta*Shp.m_GradTrial)*(ElmtSoln.m_QpGradU[1]*Shp.m_GradTest)*Ctan[0]
               +L*K*K*(Shp.m_GradTrial*Shp.m_GradTest)*Ctan[0]
               +L*d2Fdeta2*Shp.m_Trial*Shp.m_Test*Ctan[0];
    // K_eta,T
    LocalK(1,2)=L*d2FdetadT*Shp.m_Trial*Shp.m_Test*Ctan[0];
    //*************************************************************
    // K_T,eta
    LocalK(2,1)=-Latent*Shp.m_Trial*Shp.m_Test*Ctan[1];
    // K_T,T
    LocalK(2,2)=Shp.m_Trial*Shp.m_Test*Ctan[1]
               +Shp.m_GradTrial*Shp.m_GradTest*Ctan[0];

}