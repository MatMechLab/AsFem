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
//+++ Date   : 2022.10.10
//+++ Purpose: implement the residual and jacobian for stress diffusion
//+++          problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/StressDiffusionElement.h"

void StressDiffusionElement::computeAll(const FECalcType &CalcType,
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
    if(soln.m_QpU[0]||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()) {}
    if(elmtinfo.m_Dim<2){
        MessagePrinter::printErrorTxt("SmallStrainDiffusion element works only for 2d and 3d case");
        MessagePrinter::exitAsFem();
    }
    // The dofs are: 1->c, 2->ux, 3->uy, [4->uz]
    // R_c
    localR(1)=soln.m_QpV[1]*shp.m_Test
             +mate.ScalarMaterial("D")*(soln.m_QpGradU[1]*shp.m_GradTest)
             -mate.ScalarMaterial("D")*soln.m_QpU[1]*mate.ScalarMaterial("Omega")*mate.ScalarMaterial("SigmaH")*(soln.m_QpGradU[1]*shp.m_GradTest);
    // R_ux
    localR(2)=mate.Rank2Material("stress").getIthRow(1)*shp.m_GradTest;
    // R_uy
    localR(3)=mate.Rank2Material("stress").getIthRow(2)*shp.m_GradTest;
    if(elmtinfo.m_Dim==3){
        localR(4)=mate.Rank2Material("stress").getIthRow(3)*shp.m_GradTest;
    }

}
//*****************************************************************************
void StressDiffusionElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&Ctan)[3],
                                             const LocalElmtSolution &soln,
                                             const LocalShapeFun &shp,
                                             const MaterialsContainer &mate_old,
                                             const MaterialsContainer &mate,
                                             MatrixXd &LocalK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(mate_old.getScalarMaterialsNum()){}

    // K_c,c
    LocalK(1,1)=shp.m_Trial*shp.m_Test*Ctan[1]
               +mate.ScalarMaterial("D")*(shp.m_GradTrial*shp.m_GradTest)*Ctan[0]
               -mate.ScalarMaterial("D")*shp.m_Trial*mate.ScalarMaterial("Omega")*mate.ScalarMaterial("SigmaH")*(soln.m_QpGradU[1]*shp.m_GradTest)*Ctan[0]
               -mate.ScalarMaterial("D")*soln.m_QpU[1]*mate.ScalarMaterial("Omega")*mate.ScalarMaterial("SigmaH")*(shp.m_GradTrial*shp.m_GradTest)*Ctan[0];
    // K_c,ux
    LocalK(1,2)=0.0;
    // K_c,uy
    LocalK(1,3)=0.0;

    // K_ux,c
    LocalK(2,1)=mate.Rank2Material("dstressdc").getIthRow(1)*shp.m_Trial*shp.m_GradTest*Ctan[0];
    // K_ux,ux
    LocalK(2,2)=mate.Rank4Material("jacobian").getIKComponent(1,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    // K_ux,uy
    LocalK(2,3)=mate.Rank4Material("jacobian").getIKComponent(1,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

    // K_uy,c
    LocalK(3,1)=mate.Rank2Material("dstressdc").getIthRow(2)*shp.m_Trial*shp.m_GradTest*Ctan[0];
    // K_uy,ux
    LocalK(3,2)=mate.Rank4Material("jacobian").getIKComponent(2,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    // K_uy,uy
    LocalK(3,3)=mate.Rank4Material("jacobian").getIKComponent(2,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

    if(elmtinfo.m_Dim==3){
        // K_c,uz
        LocalK(1,4)=0.0;

        // K_ux,uz
        LocalK(2,4)=mate.Rank4Material("jacobian").getIKComponent(1,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

        // K_uy,uz
        LocalK(3,4)=mate.Rank4Material("jacobian").getIKComponent(2,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

        // K_uz,c
        LocalK(4,1)=mate.Rank2Material("dstressdc").getIthRow(3)*shp.m_Trial*shp.m_GradTest*Ctan[0];
        // K_uz,ux
        LocalK(4,2)=mate.Rank4Material("jacobian").getIKComponent(3,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        // K_uz,uy
        LocalK(4,3)=mate.Rank4Material("jacobian").getIKComponent(3,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        // K_uz,uz
        LocalK(4,4)=mate.Rank4Material("jacobian").getIKComponent(3,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    }

}