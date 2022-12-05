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
//+++ Date   : 2022.10.10
//+++ Purpose: implement the residual and jacobian for stress diffusion
//+++          problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/StressDiffusionElement.h"

void StressDiffusionElement::computeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::printErrorTxt("unsupported calculation type in SmallStrainDiffusionElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void StressDiffusionElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(soln.m_gpU[0]||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()) {}
    if(elmtinfo.m_dim<2){
        MessagePrinter::printErrorTxt("SmallStrainDiffusion element works only for 2d and 3d case");
        MessagePrinter::exitAsFem();
    }
    // The dofs are: 1->c, 2->ux, 3->uy, [4->uz]
    // R_c
    localR(1)=soln.m_gpV[1]*shp.m_test
             +mate.ScalarMaterial("D")*(soln.m_gpGradU[1]*shp.m_grad_test)
             -mate.ScalarMaterial("D")*soln.m_gpU[1]*mate.ScalarMaterial("Omega")*mate.ScalarMaterial("SigmaH")*(soln.m_gpGradU[1]*shp.m_grad_test);
    // R_ux
    localR(2)=mate.Rank2Material("stress").getIthRow(1)*shp.m_grad_test;
    // R_uy
    localR(3)=mate.Rank2Material("stress").getIthRow(2)*shp.m_grad_test;
    if(elmtinfo.m_dim==3){
        localR(4)=mate.Rank2Material("stress").getIthRow(3)*shp.m_grad_test;
    }

}
//*****************************************************************************
void StressDiffusionElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(mate_old.getScalarMaterialsNum()){}

    // K_c,c
    localK(1,1)=shp.m_trial*shp.m_test*ctan[1]
               +mate.ScalarMaterial("D")*(shp.m_grad_trial*shp.m_grad_test)*ctan[0]
               -mate.ScalarMaterial("D")*shp.m_trial*mate.ScalarMaterial("Omega")*mate.ScalarMaterial("SigmaH")*(soln.m_gpGradU[1]*shp.m_grad_test)*ctan[0]
               -mate.ScalarMaterial("D")*soln.m_gpU[1]*mate.ScalarMaterial("Omega")*mate.ScalarMaterial("SigmaH")*(shp.m_grad_trial*shp.m_grad_test)*ctan[0];
    // K_c,ux
    localK(1,2)=0.0;
    // K_c,uy
    localK(1,3)=0.0;

    // K_ux,c
    localK(2,1)=mate.Rank2Material("dstressdc").getIthRow(1)*shp.m_trial*shp.m_grad_test*ctan[0];
    // K_ux,ux
    localK(2,2)=mate.Rank4Material("jacobian").getIKComponent(1,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
    // K_ux,uy
    localK(2,3)=mate.Rank4Material("jacobian").getIKComponent(1,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

    // K_uy,c
    localK(3,1)=mate.Rank2Material("dstressdc").getIthRow(2)*shp.m_trial*shp.m_grad_test*ctan[0];
    // K_uy,ux
    localK(3,2)=mate.Rank4Material("jacobian").getIKComponent(2,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
    // K_uy,uy
    localK(3,3)=mate.Rank4Material("jacobian").getIKComponent(2,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

    if(elmtinfo.m_dim==3){
        // K_c,uz
        localK(1,4)=0.0;

        // K_ux,uz
        localK(2,4)=mate.Rank4Material("jacobian").getIKComponent(1,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

        // K_uy,uz
        localK(3,4)=mate.Rank4Material("jacobian").getIKComponent(2,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];

        // K_uz,c
        localK(4,1)=mate.Rank2Material("dstressdc").getIthRow(3)*shp.m_trial*shp.m_grad_test*ctan[0];
        // K_uz,ux
        localK(4,2)=mate.Rank4Material("jacobian").getIKComponent(3,1,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
        // K_uz,uy
        localK(4,3)=mate.Rank4Material("jacobian").getIKComponent(3,2,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
        // K_uz,uz
        localK(4,4)=mate.Rank4Material("jacobian").getIKComponent(3,3,shp.m_grad_test,shp.m_grad_trial)*ctan[0];
    }

}