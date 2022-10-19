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
//+++ Date   : 2021.04.10
//+++ Purpose: Calculate the double well potential 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/DoubleWellPotentialMaterial.h"

DoubleWellPotentialMaterial::DoubleWellPotentialMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
DoubleWellPotentialMaterial::~DoubleWellPotentialMaterial(){
    m_args.clean();
    m_F.clean();
    m_dFdargs.clean();
    m_d2Fdargs2.clean();
}

void DoubleWellPotentialMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void DoubleWellPotentialMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||
       mateold.getScalarMaterialsNum()){}

    //************************
    //*** here the free energy formulation is:
    //*** f=height*(x-a)^2(x-b)^2
    if(!JsonUtils::hasOnlyGivenValues(inputparams,vector<string>{"alpha","beta","w","L","eps"})){
        MessagePrinter::printErrorTxt("DoubleWellPotential material requires only: alpha, beta, w, L, and eps, "
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }
    

    m_args(1)=elmtsoln.m_gpU[1];
    computeFreeEnergyAndDerivatives(inputparams,m_args,m_F,m_dFdargs,m_d2Fdargs2);

    mate.ScalarMaterial("F")=m_F(1);
    mate.ScalarMaterial("dFdeta")=m_dFdargs(1);
    mate.ScalarMaterial("d2Fdeta2")=m_d2Fdargs2(1,1);
    mate.ScalarMaterial("eps")=JsonUtils::getValue(inputparams,"eps");
    mate.ScalarMaterial("L")=JsonUtils::getValue(inputparams,"L");
    mate.VectorMaterial("gradu")=elmtsoln.m_gpGradU[1];// the gradient of u

}
//**************************************************************
void DoubleWellPotentialMaterial::computeFreeEnergyAndDerivatives(const nlohmann::json &parameters,
                                                                  const VectorXd &args,
                                                                  VectorXd       &F,
                                                                  VectorXd       &dFdargs,
                                                                  MatrixXd       &d2Fdargs2){
    // this double well contains 1 species
    double c=args(1);
    m_a=JsonUtils::getValue(parameters,"alpha");
    m_b=JsonUtils::getValue(parameters,"beta");
    m_w=JsonUtils::getValue(parameters,"w");
    F(1)=m_w*(c-m_a)*(c-m_a)*(c-m_b)*(c-m_b);
    dFdargs(1)=2.0*m_w*(c-m_a)*(c-m_b)*(2.0*c-m_a-m_b);
    d2Fdargs2(1,1)=2.0*m_w*(m_a*m_a+4*m_a*m_b+m_b*m_b-6*m_a*c-6*m_b*c+6*c*c);

}