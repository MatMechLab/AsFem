//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Reviewer: hf @ 2023.04.04
//+++ Date    : 2022.08.20
//+++ Purpose : Implement the constitutive law for linear elastic
//+++           materials
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

    if(elmtinfo.m_dim==2){
        m_strain(3,3)=m_eps_zz;
    }

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
    double E,nu;// Youngs modulus and poisso ratio
    double K,G; // bulk moduli and shear moduli
    double lame;// lame constant

    E=nu=K=G=lame=0.0;

    if(JsonUtils::hasValue(params,"E")&&
       JsonUtils::hasValue(params,"nu")){
        E=JsonUtils::getValue(params,"E");
        nu=JsonUtils::getValue(params,"nu");
        lame=E*nu/((1.0+nu)*(1.0-2.0*nu));
        G=0.5*E/(1+nu);
    }
    else if(JsonUtils::hasValue(params,"K")&&
            JsonUtils::hasValue(params,"G")){
        K=JsonUtils::getValue(params,"K");
        G=JsonUtils::getValue(params,"G");
        E=9.0*K*G/(3.0*K+G);
        nu=(3.0*K-2.0*G)/(2.0*(3.0*K+G));
        lame=K-2.0*G/3.0;
    }
    else if(JsonUtils::hasValue(params,"Lame")&&
            JsonUtils::hasValue(params,"G")){
        lame=JsonUtils::getValue(params,"Lame");
        G=JsonUtils::getValue(params,"G");
        E=G*(3.0*lame+2.0*G)/(lame+G);
        nu=0.5*lame/(lame+G);
    }
    else{
        MessagePrinter::printErrorTxt("Invalid parameters, for linear elastic material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_sig_zz=0.0;m_eps_zz=0.0;// for plane-stress/plane-strain convertion
    if(JsonUtils::hasValue(params,"plane-strain")){
        if(dim!=2){
            MessagePrinter::printErrorTxt("plane-strain option works only for 2d case, please check your input file");
            MessagePrinter::exitAsFem();
        }
        else{
            if(JsonUtils::getBoolean(params,"plane-strain")){
                // for plane strain case
                // Within the context of tensor formulation, 
                // the condition of plane strain is inherently and automatically fulfilled.
                // Thanks hf for identifying the inconsistency present in the previous implementation..
                jacobian.setFromEAndNu(E,nu);
                stress=jacobian.doubledot(strain);
                m_eps_zz=0.0;
            }
            else{
                // for plane stress case
                // Firstly, we use the condition stress_zz=0=lame*(strain_kk)+2*mu*strain_zz to
                // obtain the correct strain_zz
                m_eps_zz=-lame*(strain(1,1)+strain(2,2))/(lame+2.0*G);
                // Next, we can get the 'correct' strain tensor from strain_zz
                m_strain_new=strain;// still 2x2 tensor
                m_strain_new(3,3)=m_eps_zz;// the correct strain tensor (3x3) for plane-stress condition
                // stress_ij=lame*delta_ij*trace(eps_ij)+2mu*eps_ij
                m_I.setToIdentity();
                m_I4Sym.setToIdentity4Symmetric();
                stress=lame*m_I*m_strain_new.trace()+2.0*G*m_strain_new;
                jacobian=m_I.otimes(m_I)*lame+m_I4Sym*2.0*G;
            }
        }
    }
    else{
        // default option is plane-strain
        if(dim==2){
            jacobian.setFromEAndNu(E,nu);
            stress=jacobian.doubledot(strain);
            m_eps_zz=0.0;
        }
        else{
            // for 3D case, just use the standard tensor formula
            jacobian.setFromEAndNu(E,nu);
            stress=jacobian.doubledot(strain);
        }
    }

}