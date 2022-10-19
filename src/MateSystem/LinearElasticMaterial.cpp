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
//+++ Date   : 2022.08.20
//+++ Purpose: Implement the constitutive law for linear elastic
//+++          materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/LinearElasticMaterial.h"

void LinearElasticMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||mate.getScalarMaterialsNum()){}
}

void LinearElasticMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    if(mateold.getRank2MaterialsNum()){}
    m_gradU.setToZeros();
    if(elmtinfo.m_dim==1){
        m_gradU.setFromGradU(elmtsoln.m_gpGradU[1]);
    }
    else if(elmtinfo.m_dim==2){
        m_gradU.setFromGradU(elmtsoln.m_gpGradU[1],elmtsoln.m_gpGradU[2]);
    }
    else if(elmtinfo.m_dim==3){
        m_gradU.setFromGradU(elmtsoln.m_gpGradU[1],elmtsoln.m_gpGradU[2],elmtsoln.m_gpGradU[3]);
    }

    computeStrain(elmtinfo.m_dim,m_gradU,m_strain);
    computeStressAndJacobian(inputparams,elmtinfo.m_dim,m_strain,m_stress,m_jacobian);

    m_devStress=m_stress.dev();
    m_devStrain=m_strain.dev();

    mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devStress.doubledot(m_devStress));
    mate.ScalarMaterial("vonMises-strain")=sqrt(1.5*m_devStrain.doubledot(m_devStrain));
    mate.ScalarMaterial("hydrostatic-stress")=m_stress.trace()/3.0;

    mate.VectorMaterial("gradux")=elmtsoln.m_gpGradU[1];
    if(elmtinfo.m_dim>=2){
        mate.VectorMaterial("graduy")=elmtsoln.m_gpGradU[2];
        if(elmtinfo.m_dim==3) mate.VectorMaterial("graduz")=elmtsoln.m_gpGradU[3];
    }
    
    mate.Rank2Material("strain")=m_strain;
    mate.Rank2Material("stress")=m_stress;
    mate.Rank4Material("jacobian")=m_jacobian;

}

void LinearElasticMaterial::computeStrain(const int &dim,const Rank2Tensor &gradU,Rank2Tensor &strain){
    if(dim){}
    strain=(gradU+gradU.transpose())*0.5;
}
void LinearElasticMaterial::computeStressAndJacobian(const nlohmann::json &params,
                                          const int &dim,
                                          const Rank2Tensor &strain,
                                          Rank2Tensor &stress,
                                          Rank4Tensor &jacobian){
    if(dim){}
    double E,nu;
    double K,G;
    double lame;

    if(JsonUtils::hasValue(params,"E")&&
       JsonUtils::hasValue(params,"nu")){
        E=JsonUtils::getValue(params,"E");
        nu=JsonUtils::getValue(params,"nu");
        jacobian.setFromEAndNu(E,nu);
    }
    else if(JsonUtils::hasValue(params,"K")&&
            JsonUtils::hasValue(params,"G")){
        K=JsonUtils::getValue(params,"K");
        G=JsonUtils::getValue(params,"G");
        jacobian.setFromKAndG(K,G);
    }
    else if(JsonUtils::hasValue(params,"Lame")&&
            JsonUtils::hasValue(params,"G")){
        lame=JsonUtils::getValue(params,"Lame");
        G=JsonUtils::getValue(params,"G");
        jacobian.setFromLameAndG(lame,G);
    }
    else{
        MessagePrinter::printErrorTxt("Invalid parameters, for linear elastic material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }

    stress=jacobian.doubledot(strain);

}