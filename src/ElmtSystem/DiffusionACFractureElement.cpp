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
//+++ Date   : 2022.11.26
//+++ Purpose: implement the residual and jacobian for diffusion 
//+++          coupled allen-cahn type phase-field fracture model
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "ElmtSystem/DiffusionACFractureElement.h"

void DiffusionACFractureElement::computeAll(const FECalcType &CalcType,
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
    else{
        MessagePrinter::printErrorTxt("unsupported calculation type in DiffusionAllenCahnFractureElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void DiffusionACFractureElement::computeResidual(const LocalElmtInfo &ElmtInfo,
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
        MessagePrinter::printErrorTxt("Diffusion coupled AllenCahnFractureElmt works only for 2d and 3d case");
        MessagePrinter::exitAsFem();
    }
    double L=Mate.ScalarMaterial("L");
    double Gc=Mate.ScalarMaterial("Gc");
    double eps=Mate.ScalarMaterial("eps");
    double Hist=Mate.ScalarMaterial("H");
    double dFdD=Mate.ScalarMaterial("dFdD");
    Rank2Tensor Stress=Mate.Rank2Material("stress");

    // R_c
    LocalR(1)=ElmtSoln.m_QpV[1]*Shp.m_Test
             +Mate.ScalarMaterial("D")*(ElmtSoln.m_QpGradU[1]*Shp.m_GradTest)
             -Mate.ScalarMaterial("D")*ElmtSoln.m_QpU[1]*Mate.ScalarMaterial("Omega")*Mate.ScalarMaterial("SigmaH")*(ElmtSoln.m_QpGradU[1]*Shp.m_GradTest);
    // R_d
    LocalR(2)=ElmtSoln.m_QpV[2]*Shp.m_Test
             +L*dg(ElmtSoln.m_QpU[2])*Hist*Shp.m_Test
             +L*dFdD*Shp.m_Test
             +L*Gc*eps*(ElmtSoln.m_QpGradU[2]*Shp.m_GradTest);
    // R_ux
    LocalR(3)=Stress.getIthRow(1)*Shp.m_GradTest;
    // R_uy
    LocalR(4)=Stress.getIthRow(2)*Shp.m_GradTest;
    if(ElmtInfo.m_Dim==3){
        // R_uz
        LocalR(5)=Stress.getIthRow(3)*Shp.m_GradTest;
    }

}
//*****************************************************************************
void DiffusionACFractureElement::computeJacobian(const LocalElmtInfo &ElmtInfo,
                                                 const double (&Ctan)[3],
                                                 const LocalElmtSolution &ElmtSoln,
                                                 const LocalShapeFun &Shp,
                                                 const MaterialsContainer &MateOld,
                                                 const MaterialsContainer &Mate,
                                                 MatrixXd &LocalK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(MateOld.getScalarMaterialsNum()){}

    int k;
    double valx,valy,valz;
    double L=Mate.ScalarMaterial("L");
    double Gc=Mate.ScalarMaterial("Gc");
    double eps=Mate.ScalarMaterial("eps");
    double Hist=Mate.ScalarMaterial("H");
    double d2FdD2=Mate.ScalarMaterial("d2FdD2");
    Rank2Tensor dStressdC=Mate.Rank2Material("dstressdc");
    Rank2Tensor dStressdD=Mate.Rank2Material("dstressdD");
    Rank2Tensor dHdstrain=Mate.Rank2Material("dHdstrain");

    // K_c,c
    LocalK(1,1)=Shp.m_Trial*Shp.m_Test*Ctan[1]
               +Mate.ScalarMaterial("D")*(Shp.m_GradTrial*Shp.m_GradTest)*Ctan[0]
               -Mate.ScalarMaterial("D")*Shp.m_Trial*Mate.ScalarMaterial("Omega")*Mate.ScalarMaterial("SigmaH")*(ElmtSoln.m_QpGradU[1]*Shp.m_GradTest)*Ctan[0]
               -Mate.ScalarMaterial("D")*ElmtSoln.m_QpU[1]*Mate.ScalarMaterial("Omega")*Mate.ScalarMaterial("SigmaH")*(Shp.m_GradTrial*Shp.m_GradTest)*Ctan[0];
    // K_c,d
    LocalK(1,2)=0.0;
    // K_c,ux
    LocalK(1,3)=0.0;
    // K_c,uy
    LocalK(1,4)=0.0;

    // K_d,c
    LocalK(2,1)=0.0;
    // K_d,d
    LocalK(2,2)=Shp.m_Trial*Shp.m_Test*Ctan[1]
               +L*d2g(ElmtSoln.m_QpU[2])*Shp.m_Trial*Hist*Shp.m_Test*Ctan[0]
               +L*d2FdD2*Shp.m_Trial*Shp.m_Test*Ctan[0]
               +L*Gc*eps*(Shp.m_GradTrial*Shp.m_GradTest)*Ctan[0];
    //*******************************
    valx=0.0;valy=0.0;valz=0.0;
    if(Mate.BooleanMaterial("finite-strain")){
        // for finite strain case
        for(k=1;k<=3;k++){
            valx+=dHdstrain(1,k)*Shp.m_GradTrial(k);
            valy+=dHdstrain(2,k)*Shp.m_GradTrial(k);
            valz+=dHdstrain(3,k)*Shp.m_GradTrial(k);
        }
    }
    else{
        // for small strain case
        for(k=1;k<=3;k++){
            valx+=0.5*(dHdstrain(1,k)+dHdstrain(k,1))*Shp.m_GradTrial(k);
            valy+=0.5*(dHdstrain(2,k)+dHdstrain(k,2))*Shp.m_GradTrial(k);
            valz+=0.5*(dHdstrain(3,k)+dHdstrain(k,3))*Shp.m_GradTrial(k);
        }
    }
    // K_d,ux
    LocalK(2,3)=L*dg(ElmtSoln.m_QpU[2])*valx*Shp.m_Test*Ctan[0];
    // K_d,uy
    LocalK(2,4)=L*dg(ElmtSoln.m_QpU[2])*valy*Shp.m_Test*Ctan[0];

    // K_ux,c
    LocalK(3,1)=dStressdC.getIthRow(1)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
    // K_ux,d
    LocalK(3,2)=dStressdD.getIthRow(1)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
    // K_ux,ux
    LocalK(3,3)=Mate.Rank4Material("jacobian").getIKComponent(1,1,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];
    // K_ux,uy
    LocalK(3,4)=Mate.Rank4Material("jacobian").getIKComponent(1,2,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

    // K_uy,c
    LocalK(4,1)=dStressdC.getIthRow(2)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
    // K_uy,d
    LocalK(4,2)=dStressdD.getIthRow(2)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
    // K_uy,ux
    LocalK(4,3)=Mate.Rank4Material("jacobian").getIKComponent(2,1,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];
    // K_uy,uy
    LocalK(4,4)=Mate.Rank4Material("jacobian").getIKComponent(2,2,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

    if(ElmtInfo.m_Dim==3){
        // K_c,uz
        LocalK(1,5)=0.0;
        
        // K_d,uz
        LocalK(2,5)=L*dg(ElmtSoln.m_QpU[2])*valz*Shp.m_Test*Ctan[0];

        // K_ux,uz
        LocalK(3,5)=Mate.Rank4Material("jacobian").getIKComponent(1,3,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

        // K_uy,uz
        LocalK(4,5)=Mate.Rank4Material("jacobian").getIKComponent(2,3,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

        // K_uz,c
        LocalK(5,1)=dStressdC.getIthRow(3)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
        // K_uz,d
        LocalK(5,2)=dStressdD.getIthRow(3)*Shp.m_Trial*Shp.m_GradTest*Ctan[0];
        // K_uz,ux
        LocalK(5,3)=Mate.Rank4Material("jacobian").getIKComponent(3,1,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];
        // K_uz,uy
        LocalK(5,4)=Mate.Rank4Material("jacobian").getIKComponent(3,2,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];
        // K_uz,uz
        LocalK(5,5)=Mate.Rank4Material("jacobian").getIKComponent(3,3,Shp.m_GradTest,Shp.m_GradTrial)*Ctan[0];

    }
}