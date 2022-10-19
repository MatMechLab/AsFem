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
//+++ Date   : 2022.10.03
//+++ Purpose: Calculate the binary mixture free energy
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BinaryMixtureMaterial.h"

BinaryMixtureMaterial::BinaryMixtureMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
BinaryMixtureMaterial::~BinaryMixtureMaterial(){
    m_args.clean();
    m_F.clean();
    m_dFdargs.clean();
    m_d2Fdargs2.clean();
}

void BinaryMixtureMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||mate.getScalarMaterialsNum()){}
}

void BinaryMixtureMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    if(elmtinfo.m_dim||mateold.getScalarMaterialsNum()) {}
    m_args(1)=elmtsoln.m_gpU[1];// 1st dof must be concentration
    computeFreeEnergyAndDerivatives(inputparams,m_args,m_F,m_dFdargs,m_d2Fdargs2);
    if(!JsonUtils::hasOnlyGivenValues(inputparams,vector<string>{"D","kappa","chi"})){
        MessagePrinter::printErrorTxt("for IdealSolutionFreeEnergyMaterial, only D, kappa, and chi are required, "
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    mate.ScalarMaterial("F")=m_F(1);
    mate.ScalarMaterial("dFdC")=m_dFdargs(1);
    mate.ScalarMaterial("d2FdC2")=m_d2Fdargs2(1,1);
    mate.ScalarMaterial("kappa")=JsonUtils::getValue(inputparams,"kappa");

    mate.ScalarMaterial("M")=JsonUtils::getValue(inputparams,"D")*m_args(1)*(1.0-m_args(1));
    mate.ScalarMaterial("dMdC")=JsonUtils::getValue(inputparams,"D")*(1.0-2.0*m_args(1));

    mate.VectorMaterial("gradc")=elmtsoln.m_gpGradU[1];// the gradient of concentration
    mate.VectorMaterial("gradmu")=elmtsoln.m_gpGradU[2];// the gradient of chemical potential
    
}
//**************************************************************
void BinaryMixtureMaterial::computeFreeEnergyAndDerivatives(const nlohmann::json &parameters,
                                                            const VectorXd &args,
                                                            VectorXd       &F,
                                                            VectorXd       &dFdargs,
                                                            MatrixXd       &d2Fdargs2){
    // binary mixture solution contains only 1 species
    double c=args(1);
    if(c<1.0e-4) c=1.0e-4;
    if(c>=1.0  ) c=1.0-1.0e-4;
    m_chi=JsonUtils::getValue(parameters,"chi");
    F(1)=c*log(c)+(1.0-c)*log(1.0-c)+m_chi*c*(1.0-c);
    dFdargs(1)=log(c)-log(1.0-c)+m_chi*(1.0-2.0*c);
    d2Fdargs2(1,1)=1.0/(1.0-c)+1.0/c-2.0*m_chi;

}