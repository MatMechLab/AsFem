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
//+++ Date   : 2022.11.19
//+++ Purpose: Calculate the material properties for diffusion coupled
//+++          elasto-plastic model in small strain case .
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/SmallStrainDiffusionJ2Material.h"

SmallStrainDiffusionJ2Material::SmallStrainDiffusionJ2Material(){
    m_args.resize(11);
}
SmallStrainDiffusionJ2Material::~SmallStrainDiffusionJ2Material(){
    m_args.clean();
}
//**********************************************************************************
void SmallStrainDiffusionJ2Material::initMaterialProperties(const nlohmann::json &inputparams,
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

void SmallStrainDiffusionJ2Material::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //*********************************************
    //*** In this model, the DoFs should be
    //*** 1.   concentration
    //*** 2.   ux
    //*** 3.   uy
    //*** [4]. uz
    //*********************************************

    mate.ScalarMaterial("D")=JsonUtils::getValue(inputparams,"D");
    mate.ScalarMaterial("Omega")=JsonUtils::getValue(inputparams,"Omega");

    if(elmtinfo.m_dim==1){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[2]);
    }
    else if(elmtinfo.m_dim==2){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[2],elmtsoln.m_gpGradU[3]);
    }
    else if(elmtinfo.m_dim==3){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[2],elmtsoln.m_gpGradU[3],elmtsoln.m_gpGradU[4]);
    }
    else{
        MessagePrinter::printErrorTxt("dim>3 is invalid for SmallStrainDiffusionJ2Material, please check your code");
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

    mate.VectorMaterial("gradc")=elmtsoln.m_gpGradU[1];

}
//***********************************************************************************
double SmallStrainDiffusionJ2Material::computeYieldFunction(const nlohmann::json &parameters,
                                        const VectorXd &args,
                                        const MaterialsContainer &mate){
    double H,YieldStress;
    YieldStress=JsonUtils::getValue(parameters,"Yield-stress");
    H=JsonUtils::getValue(parameters,"Hardening-modulus");
    m_eff_plastic_strain=args(1);// effective plastic strain
    return mate.Rank2Material("Strial").norm()-sqrt(2.0/3.0)*(YieldStress+m_eff_plastic_strain*H);
}
double SmallStrainDiffusionJ2Material::computeYieldFunctionDeriv(const nlohmann::json &parameters,
                                             const VectorXd &args,
                                             const MaterialsContainer &mate){
    if(parameters.size()||args.getM()||mate.getRank2MaterialsNum()){}
    return 0.0;
}
void SmallStrainDiffusionJ2Material::computeAdmissibleStressState(const nlohmann::json &parameters,
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
        MessagePrinter::printErrorTxt("Invalid parameters, for small strain diffusion J2 plasticity material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }
    m_I.setToIdentity();
    // for the concentration induced eigen strain
    m_Omega=JsonUtils::getValue(parameters,"Omega");
    m_Cref=JsonUtils::getValue(parameters,"Cref");

    m_I.setToIdentity();
    m_eigenstrain=m_I*(elmtsoln.m_gpU[1]-m_Cref)*m_Omega/3.0;
    m_deigenstrain_dc=m_I*(m_Omega/3.0);
    m_dtotal_strain_dc=m_deigenstrain_dc;
    // implementation for BOX 3.2(P124)
    m_plastic_strain_old=mateold.Rank2Material("plastic-strain");

    m_dev_strain=(total_strain-m_eigenstrain)-m_I*(total_strain-m_eigenstrain).trace()*(1.0/3.0);
    m_ddev_strain_dc=m_deigenstrain_dc*(-1.0)+m_I*m_deigenstrain_dc.trace()*(1.0/3.0);

    m_stress_trial=(m_dev_strain-m_plastic_strain_old)*2.0*G;
    m_dstress_trial_dc=m_ddev_strain_dc*2.0*G;

    m_args(1)=mateold.ScalarMaterial("effective-plastic-strain");// from previous step
    mate.Rank2Material("Strial")=m_stress_trial;

    m_F=computeYieldFunction(parameters,m_args,mate);
    if(m_F<=0.0){
        // for elastic case
        mate.ScalarMaterial("effective-plastic-strain")=mateold.ScalarMaterial("effective-plastic-strain");
        mate.Rank2Material("plastic-strain")=mateold.Rank2Material("plastic-strain");
        mate.Rank2Material("stress")=m_I*K*(total_strain-m_eigenstrain).trace()+m_stress_trial;
        mate.Rank2Material("dstressdc")=m_I*K*-1.0*m_dtotal_strain_dc.trace()+m_dstress_trial_dc;
        mate.Rank4Material("jacobian").setFromEAndNu(E,nu);
    }
    else{
        // for plastic case, in this linear hardening case, you don't need to do the iteration!
        // for the details, please refer to BOX 3.1 and Eq(3.3.3) and Eq(3.3.6)
        m_gamma=m_F/(2.0*G+2.0*H/3.0);
        m_N=m_stress_trial*(1.0/m_stress_trial.norm());
        m_dN_dc=m_dstress_trial_dc*(1.0/m_stress_trial.norm())
               -m_stress_trial*(m_dstress_trial_dc.doubledot(m_stress_trial)/pow(m_stress_trial.norm(),3));

        m_theta=1.0-2.0*G*m_gamma*(1.0/m_stress_trial.norm());
        m_theta_bar=1.0/(1.0+H/(3.0*G))-(1.0-m_theta);

        m_I.setToIdentity();
        m_I4Sym.setToIdentity4Symmetric();

        // update stress and strain
        mate.ScalarMaterial("effective-plastic-strain")=mateold.ScalarMaterial("effective-plastic-strain")+sqrt(2.0/3.0)*m_gamma;
        mate.Rank2Material("plastic-strain")=m_plastic_strain_old+m_N*m_gamma;
        mate.Rank2Material("stress")=m_I*(total_strain-m_eigenstrain).trace()*K
                                    +m_stress_trial-m_N*2.0*G*m_gamma;
        mate.Rank2Material("dstressdc")=m_I*m_dtotal_strain_dc.trace()*-1.0*K
                                       +m_dstress_trial_dc-m_dN_dc*2.0*G*m_gamma;
        // setup elastoplastic jacobian
        mate.Rank4Material("jacobian")=m_I.otimes(m_I)*K
                                      +(m_I4Sym-m_I.otimes(m_I)*(1.0/3.0))*2.0*G*m_theta
                                      -m_N.otimes(m_N)*2.0*G*m_theta_bar;
    }

    mate.ScalarMaterial("SigmaH")=mate.Rank2Material("dstressdc").trace()/3.0;
    mate.Rank2Material("strain")=total_strain;
    mate.Rank2Material("elastic-strain")=total_strain-mate.Rank2Material("plastic-strain")-m_eigenstrain;
}
