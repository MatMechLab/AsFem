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
//+++ Date   : 2023.07.27
//+++ Purpose: Calculate User-1 defined materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/User1Material.h"

User1Material::User1Material(){
    srand(time(00));// for random seed
}
void User1Material::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void User1Material::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]||
       mateold.getScalarMaterialsNum()||mate.getScalarMaterialsNum()){}

    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    mate.ScalarMaterial("sigma")=JsonUtils::getValue(inputparams,"sigma");// sigma
    mate.ScalarMaterial("dsigmadu")=0.0;// dsigma/dphi
    mate.ScalarMaterial("f")=JsonUtils::getValue(inputparams,"f");// F
    mate.ScalarMaterial("dfdu")=0.0;// dF/dphi
    mate.VectorMaterial("gradu")=elmtsoln.m_QpGradU[1];// the gradient of u

    mate.ScalarMaterial("myx")=1.0*rand()/RAND_MAX;
    if(elmtinfo.m_ElmtID==elmtinfo.m_ElmtsNum){
        if(elmtinfo.m_QpointID==elmtinfo.m_QpointsNum){
            mate.ScalarMaterial("myx")=1000.0;
        }
    }

}
