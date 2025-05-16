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

void MieheFractureElement::computeAll(const FECalcType &CalcType,
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
        MessagePrinter::printErrorTxt("unsupported calculation type in MieheFractureElement, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void MieheFractureElement::computeResidual(const LocalElmtInfo &ElmtInfo,
                                           const LocalElmtSolution &Soln,
                                           const LocalShapeFun &Shp,
                                           const MaterialsContainer &MateOld,
                                           const MaterialsContainer &Mate,
                                           VectorXd &LocalR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(MateOld.getScalarMaterialsNum()) {}

    if(ElmtInfo.m_Dim<2){
        MessagePrinter::printErrorTxt("MieheFractureElement works only for 2d and 3d case");
        MessagePrinter::exitAsFem();
    }
    double viscosity=Mate.ScalarMaterial("viscosity");
    double Gc=Mate.ScalarMaterial("Gc");
    double eps=Mate.ScalarMaterial("eps");
    double Hist=Mate.ScalarMaterial("H");
    double dFdD=Mate.ScalarMaterial("dFdD");
    Rank2Tensor Stress=Mate.Rank2Material("stress");
    // R_d
    LocalR(1)=viscosity*Soln.m_QpV[1]*Shp.m_Test
             +dg(Soln.m_QpU[1])*Hist*Shp.m_Test
             +dFdD*Shp.m_Test
             +Gc*eps*(Soln.m_QpGradU[1]*Shp.m_GradTest);
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
void MieheFractureElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&Ctan)[3],
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
    localK(1,1)=viscosity*shp.m_Trial*shp.m_Test*Ctan[1]
               +d2g(soln.m_QpU[1])*shp.m_Trial*Hist*shp.m_Test*Ctan[0]
               +d2FdD2*shp.m_Trial*shp.m_Test*Ctan[0]
               +Gc*eps*(shp.m_GradTrial*shp.m_GradTest)*Ctan[0];
    //*******************************
    valx=0.0;valy=0.0;valz=0.0;
    if(mate.BooleanMaterial("finite-strain")){
        // for finite strain case
        for(k=1;k<=3;k++){
            valx+=dHdstrain(1,k)*shp.m_GradTrial(k);
            valy+=dHdstrain(2,k)*shp.m_GradTrial(k);
            valz+=dHdstrain(3,k)*shp.m_GradTrial(k);
        }
    }
    else{
        // for small strain case
        for(k=1;k<=3;k++){
            valx+=0.5*(dHdstrain(1,k)+dHdstrain(k,1))*shp.m_GradTrial(k);
            valy+=0.5*(dHdstrain(2,k)+dHdstrain(k,2))*shp.m_GradTrial(k);
            valz+=0.5*(dHdstrain(3,k)+dHdstrain(k,3))*shp.m_GradTrial(k);
        }
    }
    // K_d,ux
    localK(1,2)=dg(soln.m_QpU[1])*valx*shp.m_Test*Ctan[0];
    // K_d,uy
    localK(1,3)=dg(soln.m_QpU[1])*valy*shp.m_Test*Ctan[0];

    // K_ux,d
    localK(2,1)=dStressdD.getIthRow(1)*shp.m_Trial*shp.m_GradTest*Ctan[0];
    // K_ux,ux
    localK(2,2)=mate.Rank4Material("jacobian").getIKComponent(1,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    // K_ux,uy
    localK(2,3)=mate.Rank4Material("jacobian").getIKComponent(1,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

    // K_uy,d
    localK(3,1)=dStressdD.getIthRow(2)*shp.m_Trial*shp.m_GradTest*Ctan[0];
    // K_uy,ux
    localK(3,2)=mate.Rank4Material("jacobian").getIKComponent(2,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    // K_uy,uy
    localK(3,3)=mate.Rank4Material("jacobian").getIKComponent(2,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

    if(elmtinfo.m_Dim==3){
        // K_d,uz
        localK(1,4)=dg(soln.m_QpU[1])*valz*shp.m_Test*Ctan[0];

        // K_ux,uz
        localK(2,4)=mate.Rank4Material("jacobian").getIKComponent(1,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

        // K_uy,uz
        localK(3,4)=mate.Rank4Material("jacobian").getIKComponent(2,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

        // K_uz,d
        localK(4,1)=dStressdD.getIthRow(3)*shp.m_Trial*shp.m_GradTest*Ctan[0];
        // K_uz,ux
        localK(4,2)=mate.Rank4Material("jacobian").getIKComponent(3,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        // K_uz,uy
        localK(4,3)=mate.Rank4Material("jacobian").getIKComponent(3,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        // K_uz,uz
        localK(4,4)=mate.Rank4Material("jacobian").getIKComponent(3,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

    }
}