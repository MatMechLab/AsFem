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
//+++ Date   : 2022.08.26
//+++ Purpose: Calculate the material properties required by Poisson
//+++          element. In this code, we can define:
//+++           1) Sigma
//+++           2) dSigma/du(=0)
//+++           3) F
//+++           4) dF/du(=0)
//+++          Standard benchmark test for 1d poisson equation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/Poisson1DBenchmarkMaterial.h"

void Poisson1DBenchmarkMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void Poisson1DBenchmarkMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||
       mateold.getScalarMaterialsNum()||mate.getScalarMaterialsNum()){}

    double sigma,s,a,c0,c1,x;
    
    sigma=JsonUtils::getValue(inputparams,"sigma");
    a=JsonUtils::getValue(inputparams,"a");

    c0=a/12.0-0.5;
    c1=0.0;

    x=elmtinfo.m_gpCoords0(1);
    s=1.0-a*x*x;// the source term

    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    mate.ScalarMaterial("sigma")=sigma;// sigma
    mate.ScalarMaterial("dsigmadu")=0.0;// dsigma/dphi
    mate.ScalarMaterial("f")=s;// F
    mate.ScalarMaterial("dfdu")=0.0;// dF/dphi
    mate.VectorMaterial("gradu")=elmtsoln.m_gpGradU[1];// the gradient of u
    mate.ScalarMaterial("exactsolution")=0.5*x*x-a*x*x*x*x/12.0+c0*x+c1;

}