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
//+++ Purpose: Calculate the material properties required by diffusion
//+++          element. In this code, we can define:
//+++           1) D
//+++           2) dD/dc(=0)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/ConstDiffusionMaterial.h"

void ConstDiffusionMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void ConstDiffusionMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                                       const LocalElmtInfo &elmtinfo,
                                                       const LocalElmtSolution &elmtsoln,
                                                       const MaterialsContainer &mateold,
                                                       MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]||
       mateold.getScalarMaterialsNum()||mate.getScalarMaterialsNum()){}

    mate.ScalarMaterial("D")=JsonUtils::getValue(inputparams,"D");// diffusivity
    mate.ScalarMaterial("dDdc")=0.0;// dD/dc
    mate.VectorMaterial("gradc")=elmtsoln.m_QpGradU[1];// the gradient of concentration

}