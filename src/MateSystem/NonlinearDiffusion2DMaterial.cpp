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
//+++ Purpose: Calculate the material properties required by diffusion
//+++          element. In this code, we can define:
//+++           1) D
//+++           2) dD/dc
//+++          Standard benchmark test for 2d nonlinear diffusion equation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/NonlinearDiffusion2DMaterial.h"

void NonlinearDiffusion2DMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void NonlinearDiffusion2DMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]||
       mateold.getScalarMaterialsNum()||mate.getScalarMaterialsNum()){}

    double D0,delta,x,y,c;
    
    D0=JsonUtils::getValue(inputparams,"D");
    delta=JsonUtils::getValue(inputparams,"Delta");

    x=elmtinfo.m_QpCoords0(1);
    y=elmtinfo.m_QpCoords0(2);

    c=elmtsoln.m_QpU[1];

    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    mate.ScalarMaterial("D")=D0*(1.0+c*c)*(1.0+sin(x*y)*delta);// D
    mate.ScalarMaterial("dDdc")=D0*2.0*c*(1.0+sin(x*y)*delta);// dD/dc
    mate.VectorMaterial("gradc")=elmtsoln.m_QpGradU[1];// the gradient of c

}