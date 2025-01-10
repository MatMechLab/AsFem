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
//+++ Date   : 2022.08.26
//+++ Purpose: Calculate the material properties required by Poisson
//+++          element. In this code, we can define:
//+++           1) Sigma
//+++           2) dSigma/du(=0)
//+++           3) F
//+++           4) dF/du(=0)
//+++          Standard benchmark test for 2d nonlinear poisson equation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/NonlinearPoisson3DMaterial.h"

void NonlinearPoisson3DMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void NonlinearPoisson3DMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]||
       mateold.getScalarMaterialsNum()||mate.getScalarMaterialsNum()){}

    double sigma,f,x,y,z,u;
    
    sigma=JsonUtils::getValue(inputparams,"sigma");

    x=elmtinfo.m_QpCoords0(1);
    y=elmtinfo.m_QpCoords0(2);
    z=elmtinfo.m_QpCoords0(3);

    f=sin(x*y*z);
    u=elmtsoln.m_QpU[1];

    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    mate.ScalarMaterial("sigma")=sigma*(1.2+sin(u));// sigma
    mate.ScalarMaterial("dsigmadu")=sigma*cos(u);// dsigma/dphi
    mate.ScalarMaterial("f")=f;// F
    mate.ScalarMaterial("dfdu")=0.0;// dF/dphi
    mate.VectorMaterial("gradu")=elmtsoln.m_QpGradU[1];// the gradient of u

}