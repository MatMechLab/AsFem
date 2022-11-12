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
//+++ Date   : 2022.11.12
//+++ Purpose: Calculate the residual and jacobian for dendrite growth model
//+++          based on Kobayashi's paper
//+++ Ref    : Modeling and numerical simulations of dendritic crystal growth
//+++ DOI    : https://doi.org/10.1016/0167-2789(93)90120-P
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/KobayashiElement.h"

void KobayashiElement::computeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
            const LocalElmtSolution &soln,const LocalShapeFun &shp,
            const MaterialsContainer &mate_old,const MaterialsContainer &mate,
            MatrixXd &localK,VectorXd &localR) {
    if(calctype==FECalcType::COMPUTERESIDUAL){
        computeResidual(elmtinfo,soln,shp,mate_old,mate,localR);
    }
    else if(calctype==FECalcType::COMPUTEJACOBIAN){
        computeJacobian(elmtinfo,ctan,soln,shp,mate_old,mate,localK);
    }
    else{
        MessagePrinter::printErrorTxt("unsupported calculation type in KobayashiElement, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void KobayashiElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(mate_old.getScalarMaterialsNum()) {}

    if(elmtinfo.m_dim!=2){
        MessagePrinter::printErrorTxt("Kobayashi element only works for 2d case");
        MessagePrinter::exitAsFem();
    }

    L=mate.ScalarMaterial("L");
    K=mate.ScalarMaterial("K");
    dK=mate.ScalarMaterial("dK");
    dFdeta=mate.ScalarMaterial("dFdeta");
    Latent=mate.ScalarMaterial("Latent-heat");
    // For R_eta
    V(1)=-soln.m_gpGradU[1](2);
    V(2)= soln.m_gpGradU[1](1);
    V(3)= 0.0;
    localR(1)=soln.m_gpV[1]*shp.m_test
             +L*K*dK*(V*shp.m_grad_test)
             +L*K*K*(soln.m_gpGradU[1]*shp.m_grad_test)
             +L*dFdeta*shp.m_test;
    // For R_T
    localR(2)=soln.m_gpV[2]*shp.m_test
             +soln.m_gpGradU[2]*shp.m_grad_test
             -Latent*soln.m_gpV[1]*shp.m_test;

}
//*****************************************************************************
void KobayashiElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_dt||mate_old.getScalarMaterialsNum()){}
    L=mate.ScalarMaterial("L");
    K=mate.ScalarMaterial("K");
    dK=mate.ScalarMaterial("dK");
    Latent=mate.ScalarMaterial("Latent-heat");
    d2Fdeta2=mate.ScalarMaterial("d2Fdeta2");
    d2FdetadT=mate.ScalarMaterial("d2FdetadT");
    dKdGradEta=mate.VectorMaterial("dKdGradEta");
    ddKdGradEta=mate.VectorMaterial("ddKdGradEta");
    //************************************************
    V(1)=-soln.m_gpGradU[1](2);
    V(2)= soln.m_gpGradU[1](1);
    V(3)= 0.0;
    //************************
    dV(1)=-shp.m_grad_trial(2);
    dV(2)= shp.m_grad_trial(1);
    dV(3)= 0.0;
    // K_eta,eta
    localK(1,1)=shp.m_trial*shp.m_test*ctan[1]
               +L*(dKdGradEta*shp.m_grad_trial)*dK*(V*shp.m_grad_test)*ctan[0]
               +L*K*(ddKdGradEta*shp.m_grad_trial)*(V*shp.m_grad_test)*ctan[0]
               +L*K*dK*(dV*shp.m_grad_test)*ctan[0]
               +L*2*K*(dKdGradEta*shp.m_grad_trial)*(soln.m_gpGradU[1]*shp.m_grad_test)*ctan[0]
               +L*K*K*(shp.m_grad_trial*shp.m_grad_test)*ctan[0]
               +L*d2Fdeta2*shp.m_trial*shp.m_test*ctan[0];
    // K_eta,T
    localK(1,2)=L*d2FdetadT*shp.m_trial*shp.m_test*ctan[0];
    //*************************************************************
    // K_T,eta
    localK(2,1)=-Latent*shp.m_trial*shp.m_test*ctan[1];
    // K_T,T
    localK(2,2)=shp.m_trial*shp.m_test*ctan[1]
               +shp.m_grad_trial*shp.m_grad_test*ctan[0];

}