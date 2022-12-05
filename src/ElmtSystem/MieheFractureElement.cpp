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
//+++ Date   : 2022.10.25
//+++ Purpose: implement the residual and jacobian for phase-field fracture 
//             model based on Prof. Miehe's CMAME paper
//+++          1) dD/dt=-L*(delta f/delta D)
//+++          2) div(\Sigma)=0
//+++ Ref    : A phase field model for rate-independent crack propagation: 
//             Robust algorithmic implementation based on operator splits
//+++ DOI    : https://doi.org/10.1016/j.cma.2010.04.011
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/MieheFractureElement.h"

void MieheFractureElement::computeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::printErrorTxt("unsupported calculation type in MieheFractureElement, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void MieheFractureElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(mate_old.getScalarMaterialsNum()) {}

    if(elmtinfo.m_dim<2){
        MessagePrinter::printErrorTxt("MieheFractureElement works only for 2d and 3d case");
        MessagePrinter::exitAsFem();
    }
    double viscosity=mate.ScalarMaterial("viscosity");
    double Gc=mate.ScalarMaterial("Gc");
    double eps=mate.ScalarMaterial("eps");
    double Hist=mate.ScalarMaterial("H");
    double dFdD=mate.ScalarMaterial("dFdD");
    Rank2Tensor Stress=mate.Rank2Material("stress");
    // R_d
    localR(1)=viscosity*soln.m_gpV[1]*shp.m_test
             +dg(soln.m_gpU[1])*Hist*shp.m_test
             +dFdD*shp.m_test
             +Gc*eps*(soln.m_gpGradU[1]*shp.m_grad_test);
    // R_ux
    localR(2)=Stress.getIthRow(1)*shp.m_grad_test;
    // R_uy
    localR(3)=Stress.getIthRow(2)*shp.m_grad_test;
    if(elmtinfo.m_dim==3){
        // R_uz
        localR(4)=Stress.getIthRow(3)*shp.m_grad_test;
    }

}
//*****************************************************************************
void MieheFractureElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(mate_old.getScalarMaterialsNum()){}

    int k;
    double valx,valy,valz;
    double viscosity=mate.ScalarMaterial("viscosity");
    double Gc=mate.ScalarMaterial("Gc");
    double eps=mate.ScalarMaterial("eps");
    double Hist=mate.ScalarMaterial("H");
    double d2FdD2=mate.ScalarMaterial("d2FdD2");
    Rank2Tensor dStressdD=mate.Rank2Material("dstressdD");
    Rank2Tensor dHdstrain=mate.Rank2Material("dHdstrain");

    // K_d,d
    localK(1,1)=viscosity*shp.m_trial*shp.m_test*ctan[1]
               +d2g(soln.m_gpU[1])*shp.m_trial*Hist*shp.m_test*ctan[0]
               +d2FdD2*shp.m_trial*shp.m_test*ctan[0]
               +Gc*eps*(shp.m_grad_trial*shp.m_grad_test)*ctan[0];
    //*******************************
    valx=0.0;valy=0.0;valz=0.0;
    if(mate.BooleanMaterial("finite-strain")){
        // for finite strain case
        for(k=1;k<=3;k++){
            valx+=dHdstrain(1,k)*shp.m_grad_trial(k);
            valy+=dHdstrain(2,k)*shp.m_grad_trial(k);
            valz+=dHdstrain(3,k)*shp.m_grad_trial(k);
        }
    }
    else{
        // for small strain case
        for(k=1;k<=3;k++){
            valx+=0.5*(dHdstrain(1,k)+dHdstrain(k,1))*shp.m_grad_trial(k);
            valy+=0.5*(dHdstrain(2,k)+dHdstrain(k,2))*shp.m_grad_trial(k);
            valz+=0.5*(dHdstrain(3,k)+dHdstrain(k,3))*shp.m_grad_trial(k);
        }
    }
    // K_d,ux
    localK(1,2)=dg(soln.m_gpU[1])*valx*shp.m_test*ctan[0];
    // K_d,uy
    localK(1,3)=dg(soln.m_gpU[1])*valy*shp.m_test*ctan[0];

    // K_ux,d
    localK(2,1)=dStressdD.getIthRow(1)*shp.m_trial*shp.m_grad_test*ctan[0];
    // K_ux,ux
    localK(2,2)=mate.Rank4Material("jacobian").getIKComponent(1,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
    // K_ux,uy
    localK(2,3)=mate.Rank4Material("jacobian").getIKComponent(1,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

    // K_uy,d
    localK(3,1)=dStressdD.getIthRow(2)*shp.m_trial*shp.m_grad_test*ctan[0];
    // K_uy,ux
    localK(3,2)=mate.Rank4Material("jacobian").getIKComponent(2,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
    // K_uy,uy
    localK(3,3)=mate.Rank4Material("jacobian").getIKComponent(2,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

    if(elmtinfo.m_dim==3){
        // K_d,uz
        localK(1,4)=dg(soln.m_gpU[1])*valz*shp.m_test*ctan[0];

        // K_ux,uz
        localK(2,4)=mate.Rank4Material("jacobian").getIKComponent(1,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

        // K_uy,uz
        localK(3,4)=mate.Rank4Material("jacobian").getIKComponent(2,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

        // K_uz,d
        localK(4,1)=dStressdD.getIthRow(3)*shp.m_trial*shp.m_grad_test*ctan[0];
        // K_uz,ux
        localK(4,2)=mate.Rank4Material("jacobian").getIKComponent(3,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
        // K_uz,uy
        localK(4,3)=mate.Rank4Material("jacobian").getIKComponent(3,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
        // K_uz,uz
        localK(4,4)=mate.Rank4Material("jacobian").getIKComponent(3,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

    }
}