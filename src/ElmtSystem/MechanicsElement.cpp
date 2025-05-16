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
//+++ Date   : 2022.04.10
//+++ Purpose: implement the residual and jacobian for stress
//+++          equilibrium equation
//+++          div(P)=0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/MechanicsElement.h"

void MechanicsElement::computeAll(const FECalcType &CalcType,
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
        MessagePrinter::printErrorTxt("unsupported calculation type in MechanicsElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void MechanicsElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(soln.m_QpU[0]||shp.m_Test||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()) {}

    localR(1)=mate.Rank2Material("stress").getIthRow(1)*shp.m_GradTest;
    if(elmtinfo.m_Dim>=2){
        localR(2)=mate.Rank2Material("stress").getIthRow(2)*shp.m_GradTest;
        if(elmtinfo.m_Dim==3){
            localR(3)=mate.Rank2Material("stress").getIthRow(3)*shp.m_GradTest;
        }
    }

}
//*****************************************************************************
void MechanicsElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&Ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const MaterialsContainer &mate_old,
                                 const MaterialsContainer &mate,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_Dt||Ctan[0]||soln.m_QpU[0]||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()){}
    // for K_ux,ux
    localK(1,1)=mate.Rank4Material("jacobian").getIKComponent(1,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
    if(elmtinfo.m_Dim>=2){
        // K_ux,uy
        localK(1,2)=mate.Rank4Material("jacobian").getIKComponent(1,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        
        // K_uy,ux
        localK(2,1)=mate.Rank4Material("jacobian").getIKComponent(2,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        // K_uy,uy
        localK(2,2)=mate.Rank4Material("jacobian").getIKComponent(2,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        
        if(elmtinfo.m_Dim==3){
            // K_ux,uz
            localK(1,3)=mate.Rank4Material("jacobian").getIKComponent(1,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
            
            // K_uy,uz
            localK(2,3)=mate.Rank4Material("jacobian").getIKComponent(2,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];

            // K_uz,ux
            localK(3,1)=mate.Rank4Material("jacobian").getIKComponent(3,1,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
            // K_uz,uy
            localK(3,2)=mate.Rank4Material("jacobian").getIKComponent(3,2,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
            // K_uz,uz
            localK(3,3)=mate.Rank4Material("jacobian").getIKComponent(3,3,shp.m_GradTest,shp.m_GradTrial)*Ctan[0];
        }
    }

}