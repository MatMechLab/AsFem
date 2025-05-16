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
//+++ Date   : 2021.12.25
//+++ Purpose: implement the residual and jacobian for allen-cahn 
//+++          type phase-field fracture model
//+++          1) dD/dt=-L*(delta f/delta D)
//+++          2) div(\Sigma)=0
//+++ Ref    : A continuum phase field model for fracture, by Prof. M\"uller
//+++ DOI    : https://doi.org/10.1016/j.engfracmech.2010.08.009
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/AllenCahnFractureElement.h"

void AllenCahnFractureElement::computeAll(const FECalcType &CalcType,
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
        MessagePrinter::printErrorTxt("unsupported calculation type in AllenCahnFractureElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void AllenCahnFractureElement::computeResidual(const LocalElmtInfo &ElmtInfo,
                                               const LocalElmtSolution &ElmtSoln,
                                               const LocalShapeFun &Shp,
                                               const MaterialsContainer &MateOld,
                                               const MaterialsContainer &Mate,
                                               VectorXd &LocalR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(MateOld.getScalarMaterialsNum()) {}

    if(ElmtInfo.m_Dim<2){
        MessagePrinter::printErrorTxt("AllenCahnFractureElmt only works for 2d and 3d case");
        MessagePrinter::exitAsFem();
    }
    double L=Mate.ScalarMaterial("L");
    double Gc=Mate.ScalarMaterial("Gc");
    double eps=Mate.ScalarMaterial("Eps");
    Rank2Tensor Stress=Mate.Rank2Material("Stress");
    // R_d
    LocalR(1)=ElmtSoln.m_QpV[1]*Shp.m_Test
             +L*Mate.ScalarMaterial("dFdD")*Shp.m_Test
             +L*Gc*eps*(ElmtSoln.m_QpGradU[1]*Shp.m_GradTest);
    // R_ux
    LocalR(2)=Stress.getIthRow(1)*Shp.m_GradTest;
    // R_uy
    LocalR(3)=Stress.getIthRow(2)*Shp.m_GradTest;
    if(ElmtInfo.m_Dim==3){
        // R_uz
        LocalR(4)=Stress.getIthRow(3)*Shp.m_GradTest;
    }

}
//*****************************************************************************
void AllenCahnFractureElement::computeJacobian(const LocalElmtInfo &ElmtInfo,
                                               const double (&Ctan)[3],
                                               const LocalElmtSolution &ElmtSoln,
                                               const LocalShapeFun &Shp,
                                               const MaterialsContainer &MateOld,
                                               const MaterialsContainer &Mate,
                                               MatrixXd &LocalK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(MateOld.getScalarMaterialsNum()||ElmtSoln.m_QpU[0]){}

    int k;
    double valx,valy,valz;
    double L=Mate.ScalarMaterial("L");
    double Gc=Mate.ScalarMaterial("Gc");
    double eps=Mate.ScalarMaterial("Eps");
    double d2FdD2=Mate.ScalarMaterial("d2FdD2");
    Rank2Tensor dStressdD=Mate.Rank2Material("dStressdD");

    // K_d,d
    LocalK(1,1)=Shp.m_Trial*Shp.m_Test*Ctan[1]
                   +L*d2FdD2*Shp.m_Trial*Shp.m_Test*Ctan[0]
                   +L*Gc*eps*(Shp.m_GradTrial*Shp.m_GradTest)*Ctan[0];
    //*******************************
    valx=0.0;valy=0.0;valz=0.0;
    if(Mate.BooleanMaterial("finite-strain")){
        // for finite strain case
        for(k=1;k<=3;k++){
            valx+=Mate.Rank2Material("d2FdDdStrain")(1,k)*Shp.m_GradTrial(k);
            valy+=Mate.Rank2Material("d2FdDdStrain")(2,k)*Shp.m_GradTrial(k);
            valz+=Mate.Rank2Material("d2FdDdStrain")(3,k)*Shp.m_GradTrial(k);
        }
    }
    else{
        // for small strain case
        for(k=1;k<=3;k++){
            valx+=0.5*(Mate.Rank2Material("d2FdDdStrain")(1,k)+Mate.Rank2Material("d2FdDdStrain")(k,1))*Shp.m_GradTrial(k);
            valy+=0.5*(Mate.Rank2Material("d2FdDdStrain")(2,k)+Mate.Rank2Material("d2FdDdStrain")(k,2))*Shp.m_GradTrial(k);
            valz+=0.5*(Mate.Rank2Material("d2FdDdStrain")(3,k)+Mate.Rank2Material("d2FdDdStrain")(k,3))*Shp.m_GradTrial(k);
        }
    }
    // K_d,ux
    LocalK(1,2)=L*valx*Shp.m_Test*Ctan[0];
    // K_d,uy
    LocalK(1,3)=L*valy*Shp.m_Test*Ctan[0];

    // K_ux,d
    LocalK(2,1)=dStressdD.getIthRow(1)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
    // K_ux,ux
    LocalK(2,2)=Mate.Rank4Material("Jacobian").getIKComponent(1,1,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];
    // K_ux,uy
    LocalK(2,3)=Mate.Rank4Material("Jacobian").getIKComponent(1,2,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

    // K_uy,d
    LocalK(3,1)=dStressdD.getIthRow(2)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
    // K_uy,ux
    LocalK(3,2)=Mate.Rank4Material("Jacobian").getIKComponent(2,1,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];
    // K_uy,uy
    LocalK(3,3)=Mate.Rank4Material("Jacobian").getIKComponent(2,2,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

    if(ElmtInfo.m_Dim==3){
        // K_d,uz
        LocalK(1,4)=L*valz*Shp.m_Test*Ctan[0];

        // K_ux,uz
        LocalK(2,4)=Mate.Rank4Material("Jacobian").getIKComponent(1,3,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

        // K_uy,uz
        LocalK(3,4)=Mate.Rank4Material("Jacobian").getIKComponent(2,3,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

        // K_uz,d
        LocalK(4,1)=dStressdD.getIthRow(3)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
        // K_uz,ux
        LocalK(4,2)=Mate.Rank4Material("Jacobian").getIKComponent(3,1,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];
        // K_uz,uy
        LocalK(4,3)=Mate.Rank4Material("Jacobian").getIKComponent(3,2,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];
        // K_uz,uz
        LocalK(4,4)=Mate.Rank4Material("Jacobian").getIKComponent(3,3,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

    }
}