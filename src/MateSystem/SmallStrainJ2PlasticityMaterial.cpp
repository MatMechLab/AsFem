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
//+++ Date   : 2022.11.13
//+++ Purpose: Implement the J2 plasticity material in small strain case
//+++          the algorithm is taken from Prof. Simo's book
//+++          details can be found in BOX 3.1(P122), BOX 3.2(P124)
//+++ Ref    : Computational Inelasticity, Book by J. C. Simo and Thomas J.R. Hughes
//+++ DOI    : https://link.springer.com/book/10.1007/b98904
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/SmallStrainJ2PlasticityMaterial.h"

SmallStrainJ2PlasticityMaterial::SmallStrainJ2PlasticityMaterial(){
    m_args.resize(11);
}
SmallStrainJ2PlasticityMaterial::~SmallStrainJ2PlasticityMaterial(){
    m_args.clean();
}
//**********************************************************************************
void SmallStrainJ2PlasticityMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]){}
    mate.ScalarMaterial("effective-plastic-strain")=0.0;// store the effective plastic strain
    mate.Rank2Material("plastic-strain").setToZeros();// the plastic strain tensor
}

void SmallStrainJ2PlasticityMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    if(mateold.getScalarMaterialsNum()) {}

    if(elmtinfo.m_dim==1){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[1]);
    }
    else if(elmtinfo.m_dim==2){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[1],elmtsoln.m_gpGradU[2]);
    }
    else if(elmtinfo.m_dim==3){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[1],elmtsoln.m_gpGradU[2],elmtsoln.m_gpGradU[3]);
    }
    else{
        MessagePrinter::printErrorTxt("dim>3 is invalid for SmallStrainJ2PlasticityMaterial, please check your code");
        MessagePrinter::exitAsFem();
    }

    m_total_strain=(m_GradU+m_GradU.transpose())*0.5;
    computeAdmissibleStressState(inputparams,elmtinfo,elmtsoln,mateold,m_total_strain,mate);
    
    // for postprocess
    m_devStress=mate.Rank2Material("stress").dev();
    m_dev_strain=m_total_strain.dev();
    mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devStress.doubledot(m_devStress));
    mate.ScalarMaterial("vonMises-strain")=sqrt(1.5*m_dev_strain.doubledot(m_dev_strain));

    m_dev_strain=mate.Rank2Material("plastic-strain");
    mate.ScalarMaterial("vonMises-plastic-strain")=sqrt(1.5*m_dev_strain.doubledot(m_dev_strain));

    m_dev_strain=mate.Rank2Material("elastic-strain");
    mate.ScalarMaterial("vonMises-elastic-strain")=sqrt(1.5*m_dev_strain.doubledot(m_dev_strain));

}
//***********************************************************************************
double SmallStrainJ2PlasticityMaterial::computeYieldFunction(const nlohmann::json &parameters,
                                        const VectorXd &args,
                                        const MaterialsContainer &mate){
    double H,YieldStress;
    YieldStress=JsonUtils::getValue(parameters,"Yield-stress");
    H=JsonUtils::getValue(parameters,"Hardening-modulus");
    m_eff_plastic_strain=args(1);// effective plastic strain
    return mate.Rank2Material("Strial").norm()-sqrt(2.0/3.0)*(YieldStress+m_eff_plastic_strain*H);
}
double SmallStrainJ2PlasticityMaterial::computeYieldFunctionDeriv(const nlohmann::json &parameters,
                                             const VectorXd &args,
                                             const MaterialsContainer &mate){
    if(parameters.size()||args.getM()||mate.getRank2MaterialsNum()){}
    return 0.0;
}
void SmallStrainJ2PlasticityMaterial::computeAdmissibleStressState(const nlohmann::json &parameters,
                                              const LocalElmtInfo &elmtinfo,
                                              const LocalElmtSolution &elmtsoln,
                                              const MaterialsContainer &mateold,
                                              const Rank2Tensor &total_strain,
                                              MaterialsContainer &mate){
    if(elmtinfo.m_dim||elmtsoln.m_gpU.size()){}

    double K,E,nu,G,lame;// for elastic constants
    double H;// hardening moduli

    H=JsonUtils::getValue(parameters,"Hardening-modulus");
    K=0.0;G=0.0;

    if(JsonUtils::hasValue(parameters,"E")&&
       JsonUtils::hasValue(parameters,"nu")){
        E=JsonUtils::getValue(parameters,"E");
        nu=JsonUtils::getValue(parameters,"nu");
        lame=E*nu/((1+nu)*(1-2*nu));
        K=E/(3.0*(1.0-2.0*nu));
        G=0.5*E/(1.0+nu);
    }
    else if(JsonUtils::hasValue(parameters,"K")&&
            JsonUtils::hasValue(parameters,"G")){
        K=JsonUtils::getValue(parameters,"K");
        G=JsonUtils::getValue(parameters,"G");
        lame=K-G*2.0/3.0;
        E=9.0*K*G/(3.0*K+G);
        nu=(3.0*K-2.0*G)/(2.0*(3.0*K+G));
    }
    else if(JsonUtils::hasValue(parameters,"Lame")&&
            JsonUtils::hasValue(parameters,"mu")){
        lame=JsonUtils::getValue(parameters,"Lame");
        G=JsonUtils::getValue(parameters,"mu");
        K=lame+2.0*G/3.0;
        E=9.0*K*G/(3.0*K+G);
        nu=0.5*(3.0*K-2.0*G)/(3.0*K+G);
    }
    else if(JsonUtils::hasValue(parameters,"Lame")&&
            JsonUtils::hasValue(parameters,"G")){
        lame=JsonUtils::getValue(parameters,"Lame");
        G=JsonUtils::getValue(parameters,"G");
        K=lame+2.0*G/3.0;
        E=G*(3.0*lame+2.0*G)/(lame+G);
        nu=0.5*lame/(lame+G);
    }
    else{
        MessagePrinter::printErrorTxt("Invalid parameters, for small strain J2 plasticity material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }
    // implementation for BOX 3.2(P124)
    m_I.setToIdentity();
    m_plastic_strain_old=mateold.Rank2Material("plastic-strain");
    m_dev_strain=total_strain-m_I*total_strain.trace()*(1.0/3.0);
    m_stress_trial=(m_dev_strain-m_plastic_strain_old)*2.0*G;

    m_args(1)=mateold.ScalarMaterial("effective-plastic-strain");// from previous step
    mate.Rank2Material("Strial")=m_stress_trial;

    m_F=computeYieldFunction(parameters,m_args,mate);
    if(m_F<=0.0){
        // for elastic case
        mate.ScalarMaterial("effective-plastic-strain")=mateold.ScalarMaterial("effective-plastic-strain");
        mate.Rank2Material("plastic-strain")=mateold.Rank2Material("plastic-strain");
        mate.Rank2Material("stress")=m_I*K*total_strain.trace()+m_stress_trial;
        mate.Rank4Material("jacobian").setFromEAndNu(E,nu);
    }
    else{
        // for plastic case, in this linear hardening case, you don't need to do the iteration!
        // for the details, please refer to BOX 3.1 and Eq(3.3.3) and Eq(3.3.6)
        m_gamma=m_F/(2.0*G+2.0*H/3.0);
        m_N=m_stress_trial*(1.0/m_stress_trial.norm());

        m_theta=1.0-2.0*G*m_gamma*(1.0/m_stress_trial.norm());
        m_theta_bar=1.0/(1.0+H/(3.0*G))-(1.0-m_theta);

        m_I.setToIdentity();
        m_I4Sym.setToIdentity4Symmetric();

        // update stress and strain
        mate.ScalarMaterial("effective-plastic-strain")=mateold.ScalarMaterial("effective-plastic-strain")+sqrt(2.0/3.0)*m_gamma;
        mate.Rank2Material("plastic-strain")=m_plastic_strain_old+m_N*m_gamma;
        mate.Rank2Material("stress")=m_I*total_strain.trace()*K+m_stress_trial-m_N*2.0*G*m_gamma;
        // setup elastoplastic jacobian
        mate.Rank4Material("jacobian")=m_I.otimes(m_I)*K
                                      +(m_I4Sym-m_I.otimes(m_I)*(1.0/3.0))*2.0*G*m_theta
                                      -m_N.otimes(m_N)*2.0*G*m_theta_bar;
    }
    mate.Rank2Material("strain")=total_strain;
    mate.Rank2Material("elastic-strain")=total_strain-mate.Rank2Material("plastic-strain");
}
