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
//+++ Date   : 2021.04.10
//+++ Purpose: implement the residual and jacobian for cahn-hilliard 
//+++          equation
//+++          dc/dt=div(M*grad(u))
//+++             mu=deltaF/deltaC
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/CahnHilliardElement.h"

void CahnHilliardElement::computeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::printErrorTxt("unsupported calculation type in CahnHilliardElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void CahnHilliardElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_dt||soln.m_gpU[0]||shp.m_test||mate_old.getScalarMaterialsNum()) {}
    // R_c
    localR(1)=soln.m_gpV[1]*shp.m_test
             +mate.ScalarMaterial("M")*(soln.m_gpGradU[2]*shp.m_grad_test);
    // R_mu
    localR(2)=soln.m_gpU[2]*shp.m_test
             -mate.ScalarMaterial("dFdC")*shp.m_test
             -mate.ScalarMaterial("kappa")*soln.m_gpGradU[1]*shp.m_grad_test;

}
//*****************************************************************************
void CahnHilliardElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_dt||soln.m_gpU[0]||mate_old.getScalarMaterialsNum()){}
    // K_c,c
    localK(1,1)=shp.m_trial*shp.m_test*ctan[1]
               +mate.ScalarMaterial("dMdC")*shp.m_trial*(soln.m_gpGradU[2]*shp.m_grad_test)*ctan[0];
    // K_c,mu
    localK(1,2)=mate.ScalarMaterial("M")*(shp.m_grad_trial*shp.m_grad_test)*ctan[0];

    // K_mu,c
    localK(2,1)=-mate.ScalarMaterial("d2FdC2")*shp.m_trial*shp.m_test*ctan[0]
                -mate.ScalarMaterial("kappa")*shp.m_grad_trial*shp.m_grad_test*ctan[0];
    // K_mu,mu
    localK(2,2)=shp.m_trial*shp.m_test*ctan[0];

}