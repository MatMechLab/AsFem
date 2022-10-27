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
//+++ Date   : 2022.10.10
//+++ Purpose: Calculate the material properties required by stress-diffusion
//+++          element. In this code, we can define:
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/SmallStrainDiffusionMaterial.h"

void SmallStrainDiffusionMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void SmallStrainDiffusionMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(mateold.getScalarMaterialsNum()){}

    if(JsonUtils::hasValue(inputparams,"finite-strain")){
        mate.BooleanMaterial("finite-strain")=JsonUtils::getBoolean(inputparams,"finite-strain");
    }
    else{
        mate.BooleanMaterial("finite-strain")=false;// use small strain deformation as the default option
    }

    m_Omega=JsonUtils::getValue(inputparams,"Omega");
    m_cref=JsonUtils::getValue(inputparams,"cref");

    mate.ScalarMaterial("D")=JsonUtils::getValue(inputparams,"D");// diffusivity
    mate.ScalarMaterial("Omega")=m_Omega;// partial molar volume
    mate.VectorMaterial("gradc")=elmtsoln.m_gpGradU[1];// the gradient of concentration

    if(elmtinfo.m_dim==2){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[2],elmtsoln.m_gpGradU[3]);// grad(ux), grad(uy)
    }
    else if(elmtinfo.m_dim==3){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[2],elmtsoln.m_gpGradU[3],elmtsoln.m_gpGradU[4]);// grad(ux), grad(uy)
    }
    else{
        MessagePrinter::printErrorTxt("SmallStrainDiffusionMaterial works only for 2d and 3d case, please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_I.setToIdentity();
    computeStrain(elmtinfo.m_dim,m_GradU,m_totalstrain);
    
    m_c=elmtsoln.m_gpU[1];
    m_mechstrain=m_totalstrain-m_I*(m_c-m_cref)*m_Omega/3.0;
    m_dmechstrain_dc=m_I*(-1.0)*m_Omega/3.0;

    computeStressAndJacobian(inputparams,elmtinfo.m_dim,m_mechstrain,m_stress,m_jacobian);

    m_devstress=m_stress.dev();
    m_dstress_dc=m_jacobian.doubledot(m_dmechstrain_dc);

    mate.ScalarMaterial("D")=JsonUtils::getValue(inputparams,"D");
    mate.ScalarMaterial("Omega")=JsonUtils::getValue(inputparams,"Omega");
    mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devstress.doubledot(m_devstress));
    mate.ScalarMaterial("hydrostatic-stress")=m_stress.trace()/3.0;

    mate.ScalarMaterial("SigmaH")=m_dstress_dc.trace()/3.0;// this is the factor of hydrostatic stress, not itself !!!

    mate.Rank2Material("strain")=m_totalstrain;
    mate.Rank2Material("stress")=m_stress;
    mate.Rank2Material("dstressdc")=m_dstress_dc;

    mate.Rank4Material("jacobian")=m_jacobian;

}
//**************************************************************************
void SmallStrainDiffusionMaterial::computeStrain(const int &dim,const Rank2Tensor &gradU,Rank2Tensor &strain){
    if(dim){}
    m_I.setToIdentity();
    strain=(gradU+gradU.transpose())*0.5;// here the strain is small strain
}
void SmallStrainDiffusionMaterial::computeStressAndJacobian(const nlohmann::json &params,
                                          const int &dim,
                                          const Rank2Tensor &strain,
                                          Rank2Tensor &stress,
                                          Rank4Tensor &jacobian){
    if(dim){}
    double E=0.0,nu=0.0;
    double K=0.0,G=0.0;
    double lame=0.0;

    jacobian.setToZeros();

    if(JsonUtils::hasValue(params,"E")&&
       JsonUtils::hasValue(params,"nu")){
        E=JsonUtils::getValue(params,"E");
        nu=JsonUtils::getValue(params,"nu");
        lame=E*nu/((1+nu)*(1-2*nu));
        G=0.5*E/(1.0+nu);
        jacobian.setFromEAndNu(E,nu);
    }
    else if(JsonUtils::hasValue(params,"K")&&
            JsonUtils::hasValue(params,"G")){
        K=JsonUtils::getValue(params,"K");
        G=JsonUtils::getValue(params,"G");
        lame=K-G*2.0/3.0;
        jacobian.setFromKAndG(K,G);
    }
    else if(JsonUtils::hasValue(params,"Lame")&&
            JsonUtils::hasValue(params,"mu")){
        lame=JsonUtils::getValue(params,"Lame");
        G=JsonUtils::getValue(params,"mu");
        jacobian.setFromLameAndG(lame,G);
    }
    else{
        MessagePrinter::printErrorTxt("Invalid parameters, for small strain diffusion material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }

    stress=jacobian.doubledot(strain);

}