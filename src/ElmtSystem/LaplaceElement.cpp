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
//+++ Date   : 2022.11.20
//+++ Purpose: implement the residual and jacobian for Laplace operator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/LaplaceElement.h"

void LaplaceElement::computeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::printErrorTxt("unsupported calculation type in PoissonElmt, please check your related code");
        MessagePrinter::exitAsFem();
    }
}
//***************************************************************************
void LaplaceElement::computeResidual(const LocalElmtInfo &elmtinfo,
                                     const LocalElmtSolution &soln,
                                     const LocalShapeFun &shp,
                                     const MaterialsContainer &mate_old,
                                     const MaterialsContainer &mate,
                                     VectorXd &LocalR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.m_Dt||soln.m_QpU[0]||shp.m_Test||mate_old.getScalarMaterialsNum()||mate.getScalarMaterialsNum()) {}

    LocalR(1)=mate.ScalarMaterial("sigma")*(soln.m_QpGradU[1]*shp.m_GradTest);

}
//*****************************************************************************
void LaplaceElement::computeJacobian(const LocalElmtInfo &elmtinfo,const double (&Ctan)[3],
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
               +mate.ScalarMaterial("sigma")*(shp.m_GradTrial*shp.m_GradTest)*Ctan[0];

}