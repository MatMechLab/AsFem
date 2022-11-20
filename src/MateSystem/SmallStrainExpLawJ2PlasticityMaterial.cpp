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
//+++ Date   : 2022.11.20
//+++ Purpose: Implement the J2 plasticity material in small strain case
//+++          with exponential law hardening.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/SmallStrainExpLawJ2PlasticityMaterial.h"

SmallStrainExpLawJ2PlasticityMaterial::SmallStrainExpLawJ2PlasticityMaterial(){
    m_args.resize(11);
}
SmallStrainExpLawJ2PlasticityMaterial::~SmallStrainExpLawJ2PlasticityMaterial(){
    m_args.clean();
}
//**********************************************************************************
void SmallStrainExpLawJ2PlasticityMaterial::initMaterialProperties(const nlohmann::json &inputparams,
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

void SmallStrainExpLawJ2PlasticityMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
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
        MessagePrinter::printErrorTxt("dim>3 is invalid for SmallStrainExpLawJ2PlasticityMaterial, please check your code");
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
double SmallStrainExpLawJ2PlasticityMaterial::computeYieldFunction(const nlohmann::json &parameters,
                                        const VectorXd &args,
                                        const MaterialsContainer &mate){
    double Kinf,K0,delta,YieldStress,G,gamma;
    YieldStress=JsonUtils::getValue(parameters,"Yield-stress");
    Kinf=JsonUtils::getValue(parameters,"Kinf");
    K0=JsonUtils::getValue(parameters,"K0");
    delta=JsonUtils::getValue(parameters,"delta");

    m_eff_plastic_strain=args(1);// effective plastic strain from previous step(old)
    G=args(2);// shear moduli
    gamma=args(3);// plastic multiplier

    m_Hderiv=(Kinf-K0)*exp(-delta*m_eff_plastic_strain-delta*sqrt(2.0/3.0)*gamma)*delta;

    return mate.Rank2Material("Strial").norm()
          -2.0*G*gamma
          -sqrt(2.0/3.0)*(YieldStress
                         +(Kinf-K0)*(1.0-exp(-delta*m_eff_plastic_strain-delta*sqrt(2.0/3.0)*gamma)));
}
double SmallStrainExpLawJ2PlasticityMaterial::computeYieldFunctionDeriv(const nlohmann::json &parameters,
                                             const VectorXd &args,
                                             const MaterialsContainer &mate){
    if(mate.getScalarMaterialsNum()){}
    double Kinf,K0,delta,G,gamma;
    Kinf=JsonUtils::getValue(parameters,"Kinf");
    K0=JsonUtils::getValue(parameters,"K0");
    delta=JsonUtils::getValue(parameters,"delta");

    m_eff_plastic_strain=args(1);// effective plastic strain from previous step(old)
    G=args(2);// shear moduli
    gamma=args(3);// plastic multiplier

    m_Hderiv=(Kinf-K0)*exp(-delta*m_eff_plastic_strain-delta*sqrt(2.0/3.0)*gamma)*delta;

    return -2.0*G
           -(Kinf-K0)*exp(-delta*m_eff_plastic_strain-delta*sqrt(2.0/3.0)*gamma)*(2.0/3.0)*delta;
}
void SmallStrainExpLawJ2PlasticityMaterial::computeAdmissibleStressState(const nlohmann::json &parameters,
                                              const LocalElmtInfo &elmtinfo,
                                              const LocalElmtSolution &elmtsoln,
                                              const MaterialsContainer &mateold,
                                              const Rank2Tensor &total_strain,
                                              MaterialsContainer &mate){
    if(elmtinfo.m_dim||elmtsoln.m_gpU.size()){}

    double K,E,nu,G,lame;// for elastic constants
    double H;// hardening moduli

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
        MessagePrinter::printErrorTxt("Invalid parameters, for SmallStrainExpLawJ2PlasticityMaterial, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }
    // implementation for BOX 3.2(P124)
    m_I.setToIdentity();
    m_plastic_strain_old=mateold.Rank2Material("plastic-strain");
    m_dev_strain=total_strain-m_I*total_strain.trace()*(1.0/3.0);
    m_stress_trial=(m_dev_strain-m_plastic_strain_old)*2.0*G;

    m_args(1)=mateold.ScalarMaterial("effective-plastic-strain");// from previous step(old)
    m_args(2)=G;
    mate.Rank2Material("Strial")=m_stress_trial;

    m_args(3)=0.0;// for trial state, the gamma should be zero !
    m_F=computeYieldFunction(parameters,m_args,mate);
    if(m_F<=0.0){
        // for elastic case
        mate.ScalarMaterial("effective-plastic-strain")=mateold.ScalarMaterial("effective-plastic-strain");
        mate.Rank2Material("plastic-strain")=mateold.Rank2Material("plastic-strain");
        mate.Rank2Material("stress")=m_I*K*total_strain.trace()+m_stress_trial;
        mate.Rank4Material("jacobian").setFromEAndNu(E,nu);
    }
    else{
        // for plastic case, we do the NR iteration and find out the correct plastic multiplier
        m_N=m_stress_trial*(1.0/m_stress_trial.norm());
        m_maxiters=500;
        if(JsonUtils::hasValue(parameters,"maxiters")){
            m_maxiters=JsonUtils::getInteger(parameters,"maxiters");
            if(m_maxiters<1){
                MessagePrinter::printErrorTxt("Invalid maxiters for plastic deformation iteration, please check your input file");
                MessagePrinter::exitAsFem();
            }
        }
        m_tolerance=1.0e-4;
        if(JsonUtils::hasValue(parameters,"tolerance")){
            m_tolerance=JsonUtils::getValue(parameters,"tolerance");
            if(m_tolerance<1.0e-9){
                MessagePrinter::printErrorTxt("The tolerance(="+to_string(m_tolerance)+") is too small for the yield function in plastic deformation iteration, please increase it to a larger value");
                MessagePrinter::exitAsFem();
            }
        }
        m_iters=0;m_gamma=0.0;
        while (m_iters<m_maxiters){
            m_args(3)=m_gamma;
            m_F=computeYieldFunction(parameters,m_args,mate);
            m_dF=computeYieldFunctionDeriv(parameters,m_args,mate);
            if(m_F<m_tolerance) break;
            m_gamma+=-m_F/m_dF;
            m_iters+=1;
        }
        if(m_F>m_tolerance){
            MessagePrinter::printErrorTxt("Maximum iteration reaches for plastic deformation but |F|<tol failed, please either increase your maxiters/tolerance or check your code");
            MessagePrinter::exitAsFem();
        }
        m_args(3)=m_gamma;
        m_F=computeYieldFunction(parameters,m_args,mate);// update dH/dgamma

        H=m_Hderiv;// dH/dgamma
        m_theta=1.0/(1.0+H/(3.0*G));

        m_I.setToIdentity();
        m_I4Sym.setToIdentity4Symmetric();

        // update stress and strain
        mate.ScalarMaterial("effective-plastic-strain")=mateold.ScalarMaterial("effective-plastic-strain")+sqrt(2.0/3.0)*m_gamma;
        mate.Rank2Material("plastic-strain")=m_plastic_strain_old+m_N*m_gamma;
        mate.Rank2Material("stress")=m_I*total_strain.trace()*K+m_stress_trial-m_N*2.0*G*m_gamma;
        // setup elastoplastic jacobian
        mate.Rank4Material("jacobian")=m_I.otimes(m_I)*K
                                      +(m_I4Sym-m_I.otimes(m_I)*(1.0/3.0))*2.0*G
                                      -m_N.otimes(m_N)*2.0*G*m_theta;
    }
    mate.Rank2Material("strain")=total_strain;
    mate.Rank2Material("elastic-strain")=total_strain-mate.Rank2Material("plastic-strain");
}
