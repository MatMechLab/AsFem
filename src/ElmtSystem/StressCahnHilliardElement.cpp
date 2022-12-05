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
//+++ Date   : 2022.10.24
//+++ Purpose: implement the residual and jacobian for mechanically
//+++          coupled cahn-hilliard equation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/StressCahnHilliardElement.h"

void StressCahnHilliardElement::computeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::printErrorTxt("unsupported calculation type in StressCahnHilliardElement, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void StressCahnHilliardElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(soln.m_gpU[0]||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()) {}
    // The dofs are: 1->c, 2->mu, 3->ux, 4->uy, [5->uz]

    if(elmtinfo.m_dim<2){
        MessagePrinter::printErrorTxt("StressCahnHilliard element works only for 2d and 3d case");
        MessagePrinter::exitAsFem();
    }

    // R_c
    localR(1)=soln.m_gpV[1]*shp.m_test
             +mate.ScalarMaterial("M")*(soln.m_gpGradU[2]*shp.m_grad_test);
    // R_mu
    localR(2)=soln.m_gpU[2]*shp.m_test
             -mate.ScalarMaterial("dFdC")*shp.m_test
             -mate.ScalarMaterial("kappa")*soln.m_gpGradU[1]*shp.m_grad_test;

    // R_ux
    localR(3)=mate.Rank2Material("stress").getIthRow(1)*shp.m_grad_test;
    // R_uy
    localR(4)=mate.Rank2Material("stress").getIthRow(2)*shp.m_grad_test;
    // R_uz
    if(elmtinfo.m_dim==3) localR(5)=mate.Rank2Material("stress").getIthRow(3)*shp.m_grad_test;

}
//*****************************************************************************
void StressCahnHilliardElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(mate_old.getScalarMaterialsNum()){}

    double valx,valy,valz;
    Rank2Tensor d2FdCdStrain=mate.Rank2Material("d2FdCdStrain");
    Rank2Tensor dStressdC=mate.Rank2Material("dstressdc");

    // K_c,c
    localK(1,1)=shp.m_trial*shp.m_test*ctan[1]
               +mate.ScalarMaterial("dMdC")*shp.m_trial*(soln.m_gpGradU[2]*shp.m_grad_test)*ctan[0];
    // K_c,mu
    localK(1,2)=mate.ScalarMaterial("M")*(shp.m_grad_trial*shp.m_grad_test)*ctan[0];
    // K_c,ux
    localK(1,3)=0.0;
    // K_c,uy
    localK(1,4)=0.0;

    valx=0.0;valy=0.0;valz=0.0;
    for(int k=1;k<=3;k++){
        valx+=0.5*(d2FdCdStrain(1,k)+d2FdCdStrain(k,1))*shp.m_grad_trial(k);
        valy+=0.5*(d2FdCdStrain(2,k)+d2FdCdStrain(k,2))*shp.m_grad_trial(k);
        valz+=0.5*(d2FdCdStrain(3,k)+d2FdCdStrain(k,3))*shp.m_grad_trial(k);
    }

    // K_mu,c
    localK(2,1)=-mate.ScalarMaterial("d2FdC2")*shp.m_trial*shp.m_test*ctan[0]
                -mate.ScalarMaterial("kappa")*shp.m_grad_trial*shp.m_grad_test*ctan[0];
    // K_mu,mu
    localK(2,2)=shp.m_trial*shp.m_test*ctan[0];
    // K_mu,ux
    localK(2,3)=-1.0*valx*shp.m_test*ctan[0];
    // K_mu,uy
    localK(2,4)=-1.0*valy*shp.m_test*ctan[0];

    // K_ux,c
    localK(3,1)=dStressdC.getIthRow(1)*shp.m_trial*shp.m_grad_test*ctan[0];
    // K_ux,mu
    localK(3,2)=0.0;
    // K_ux,ux
    localK(3,3)=mate.Rank4Material("jacobian").getIKComponent(1,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
    // K_ux,uy
    localK(3,4)=mate.Rank4Material("jacobian").getIKComponent(1,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

    // K_uy,c
    localK(4,1)=dStressdC.getIthRow(2)*shp.m_trial*shp.m_grad_test*ctan[0];
    // K_uy,mu
    localK(4,2)=0.0;
    // K_uy,ux
    localK(4,3)=mate.Rank4Material("jacobian").getIKComponent(2,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
    // K_uy,uy
    localK(4,4)=mate.Rank4Material("jacobian").getIKComponent(2,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

    if(elmtinfo.m_dim==3){
        // K_c,uz
        localK(1,5)=0.0;

        // K_mu,uz
        localK(2,5)=-1.0*valz*shp.m_test*ctan[0];

        // K_ux,uz
        localK(3,5)=mate.Rank4Material("jacobian").getIKComponent(1,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

        // K_uy,uz
        localK(4,5)=mate.Rank4Material("jacobian").getIKComponent(2,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

        // K_uz,c
        localK(5,1)=dStressdC.getIthRow(3)*shp.m_trial*shp.m_grad_test*ctan[0];
        // K_uz,mu
        localK(5,2)=0.0;
        // K_uz,ux
        localK(5,3)=mate.Rank4Material("jacobian").getIKComponent(3,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
        // K_uz,uy
        localK(5,4)=mate.Rank4Material("jacobian").getIKComponent(3,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
        // K_uz,uz
        localK(5,5)=mate.Rank4Material("jacobian").getIKComponent(3,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
    }

}