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
                                 VectorXd &LocalR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(soln.m_QpU[0]||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()) {}
    // The dofs are: 1->c, 2->mu, 3->ux, 4->uy, [5->uz]

    if(elmtinfo.m_Dim<2){
        MessagePrinter::printErrorTxt("StressCahnHilliard element works only for 2d and 3d case");
        MessagePrinter::exitAsFem();
    }

    // R_c
    LocalR(1)=soln.m_QpV[1]*shp.m_Test
             +mate.ScalarMaterial("M")*(soln.m_QpGradU[2]*shp.m_GradTest);
    // R_mu
    LocalR(2)=soln.m_QpU[2]*shp.m_Test
             -mate.ScalarMaterial("dFdC")*shp.m_Test
             -mate.ScalarMaterial("kappa")*soln.m_QpGradU[1]*shp.m_GradTest;

    // R_ux
    LocalR(3)=mate.Rank2Material("stress").getIthRow(1)*shp.m_GradTest;
    // R_uy
    LocalR(4)=mate.Rank2Material("stress").getIthRow(2)*shp.m_GradTest;
    // R_uz
    if(elmtinfo.m_Dim==3) LocalR(5)=mate.Rank2Material("stress").getIthRow(3)*shp.m_GradTest;

}
//*****************************************************************************
void StressCahnHilliardElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&Ctan)[3],
                                                const LocalElmtSolution &soln,
                                                const LocalShapeFun &shp,
                                                const MaterialsContainer &mate_old,
                                                const MaterialsContainer &mate,
                                                MatrixXd &LocalK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(mate_old.getScalarMaterialsNum()){}

    double valx,valy,valz;
    Rank2Tensor d2FdCdStrain=mate.Rank2Material("d2FdCdStrain");
    Rank2Tensor dStressdC=mate.Rank2Material("dstressdc");

    // K_c,c
    LocalK(1,1)=shp.m_Trial*shp.m_Test*Ctan[1]
               +mate.ScalarMaterial("dMdC")*shp.m_Trial*(soln.m_QpGradU[2]*shp.m_GradTest)*Ctan[0];
    // K_c,mu
    LocalK(1,2)=mate.ScalarMaterial("M")*(shp.m_GradTrial*shp.m_GradTest)*Ctan[0];
    // K_c,ux
    LocalK(1,3)=0.0;
    // K_c,uy
    LocalK(1,4)=0.0;

    valx=0.0;valy=0.0;valz=0.0;
    for(int k=1;k<=3;k++){
        valx+=0.5*(d2FdCdStrain(1,k)+d2FdCdStrain(k,1))*shp.m_GradTrial(k);
        valy+=0.5*(d2FdCdStrain(2,k)+d2FdCdStrain(k,2))*shp.m_GradTrial(k);
        valz+=0.5*(d2FdCdStrain(3,k)+d2FdCdStrain(k,3))*shp.m_GradTrial(k);
    }

    // K_mu,c
    LocalK(2,1)=-mate.ScalarMaterial("d2FdC2")*shp.m_Trial*shp.m_Test*Ctan[0]
                -mate.ScalarMaterial("kappa")*shp.m_GradTrial*shp.m_GradTest*Ctan[0];
    // K_mu,mu
    LocalK(2,2)=shp.m_Trial*shp.m_Test*Ctan[0];
    // K_mu,ux
    LocalK(2,3)=-1.0*valx*shp.m_Test*Ctan[0];
    // K_mu,uy
    LocalK(2,4)=-1.0*valy*shp.m_Test*Ctan[0];

    // K_ux,c
    LocalK(3,1)=dStressdC.getIthRow(1)*shp.m_Trial*shp.m_GradTest*Ctan[0];
    // K_ux,mu
    LocalK(3,2)=0.0;
    // K_ux,ux
    LocalK(3,3)=mate.Rank4Material("jacobian").getIKComponent(1,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    // K_ux,uy
    LocalK(3,4)=mate.Rank4Material("jacobian").getIKComponent(1,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

    // K_uy,c
    LocalK(4,1)=dStressdC.getIthRow(2)*shp.m_Trial*shp.m_GradTest*Ctan[0];
    // K_uy,mu
    LocalK(4,2)=0.0;
    // K_uy,ux
    LocalK(4,3)=mate.Rank4Material("jacobian").getIKComponent(2,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    // K_uy,uy
    LocalK(4,4)=mate.Rank4Material("jacobian").getIKComponent(2,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

    if(elmtinfo.m_Dim==3){
        // K_c,uz
        LocalK(1,5)=0.0;

        // K_mu,uz
        LocalK(2,5)=-1.0*valz*shp.m_Test*Ctan[0];

        // K_ux,uz
        LocalK(3,5)=mate.Rank4Material("jacobian").getIKComponent(1,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

        // K_uy,uz
        LocalK(4,5)=mate.Rank4Material("jacobian").getIKComponent(2,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

        // K_uz,c
        LocalK(5,1)=dStressdC.getIthRow(3)*shp.m_Trial*shp.m_GradTest*Ctan[0];
        // K_uz,mu
        LocalK(5,2)=0.0;
        // K_uz,ux
        LocalK(5,3)=mate.Rank4Material("jacobian").getIKComponent(3,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        // K_uz,uy
        LocalK(5,4)=mate.Rank4Material("jacobian").getIKComponent(3,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        // K_uz,uz
        LocalK(5,5)=mate.Rank4Material("jacobian").getIKComponent(3,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    }

}