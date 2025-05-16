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
//+++ Date   : 2021.04.10
//+++ Purpose: implement the residual and jacobian for Poisson equation
//+++          Sigma*div(grad(phi))=F
//+++          Both Sigma and F could be nonlinear function !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/PoissonElement.h"

void PoissonElement::computeAll(const FECalcType &CalcType,
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
        MessagePrinter::printErrorTxt("unsupported calculation type in PoissonElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void PoissonElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                     const LocalElmtSolution &soln,
                                     const LocalShapeFun &shp,
                                     const MaterialsContainer &mate_old,
                                     const MaterialsContainer &mate,
                                     VectorXd &LocalR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_Dt||soln.m_QpU[0]||shp.m_Test||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()) {}

    LocalR(1)=mate.ScalarMaterial("sigma")*(soln.m_QpGradU[1]*shp.m_GradTest)
             +mate.ScalarMaterial("f")*shp.m_Test;

}
//*****************************************************************************
void PoissonElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&Ctan)[3],
                                     const LocalElmtSolution &soln,
                                     const LocalShapeFun &shp,
                                     const MaterialsContainer &mate_old,
                                     const MaterialsContainer &mate,
                                     MatrixXd &LocalK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_Dt||Ctan[0]||soln.m_QpU[0]||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()){}

    LocalK(1,1)=mate.ScalarMaterial("dsigmadu")*shp.m_Trial*(soln.m_QpGradU[1]*shp.m_GradTest)*Ctan[0]
               +mate.ScalarMaterial("sigma")*(shp.m_GradTrial*shp.m_GradTest)*Ctan[0]
               +mate.ScalarMaterial("dfdu")*shp.m_Trial*shp.m_Test*Ctan[0];

}