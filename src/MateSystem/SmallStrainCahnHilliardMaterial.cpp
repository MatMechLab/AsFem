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
//+++ Date   : 2022.10.24
//+++ Purpose: Calculate the material properties required by mech-cahnhilliard
//+++          element.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/SmallStrainCahnHilliardMaterial.h"

SmallStrainCahnHilliardMaterial::SmallStrainCahnHilliardMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
SmallStrainCahnHilliardMaterial::~SmallStrainCahnHilliardMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
//*******************************************************
void SmallStrainCahnHilliardMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void SmallStrainCahnHilliardMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(mateold.getScalarMaterialsNum()){}

    m_Omega=JsonUtils::getValue(inputparams,"Omega");
    m_cref=JsonUtils::getValue(inputparams,"cref");

    mate.ScalarMaterial("Omega")=m_Omega;// partial molar volume
    mate.ScalarMaterial("kappa")=JsonUtils::getValue(inputparams,"kappa");;
    mate.VectorMaterial("gradc")=elmtsoln.m_gpGradU[1];// the gradient of concentration

    mate.ScalarMaterial("M")=JsonUtils::getValue(inputparams,"D");// for mobility
    mate.ScalarMaterial("dMdC")=0.0;

    if(elmtinfo.m_dim==2){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[3],elmtsoln.m_gpGradU[4]);// grad(ux), grad(uy)
    }
    else if(elmtinfo.m_dim==3){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[3],elmtsoln.m_gpGradU[4],elmtsoln.m_gpGradU[5]);// grad(ux), grad(uy)
    }
    else{
        MessagePrinter::printErrorTxt("SmallStrainCahnHilliardMaterial works only for 2d and 3d case, please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_I.setToIdentity();
    computeStrain(elmtinfo.m_dim,m_GradU,m_totalstrain);
    
    m_c=elmtsoln.m_gpU[1];// concentration
    m_mechstrain=m_totalstrain-m_I*(m_c-m_cref)*m_Omega/3.0;
    m_dmechstrain_dc=m_I*(-1.0)*m_Omega/3.0;

    computeStressAndJacobian(inputparams,elmtinfo.m_dim,m_mechstrain,m_stress,m_jacobian);

    m_devstress=m_stress.dev();
    m_dstress_dc=m_jacobian.doubledot(m_dmechstrain_dc);

    m_args(1)=m_c;// concentration
    computeFreeEnergyAndDerivatives(inputparams,m_args,m_F,m_dFdargs,m_d2Fdargs2);

    // for free energies
    mate.ScalarMaterial("F")=m_F(1)+m_stress.doubledot(m_mechstrain)*0.5;
    mate.ScalarMaterial("dFdC")=m_dFdargs(1)+m_stress.doubledot(m_dmechstrain_dc);
    mate.ScalarMaterial("d2FdC2")=m_d2Fdargs2(1,1)+m_dstress_dc.doubledot(m_dmechstrain_dc);
    mate.Rank2Material("d2FdCdStrain")=m_dstress_dc;

    mate.Rank2Material("strain")=m_totalstrain;
    mate.Rank2Material("stress")=m_stress;
    mate.Rank2Material("dstressdc")=m_dstress_dc;

    mate.Rank4Material("jacobian")=m_jacobian;

    // for auxilary variables
    mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devstress.doubledot(m_devstress));
    mate.ScalarMaterial("hydrostatic-stress")=m_stress.trace()/3.0;

}
//**************************************************************************
void SmallStrainCahnHilliardMaterial::computeStrain(const int &dim,const Rank2Tensor &gradU,Rank2Tensor &strain){
    if(dim){}
    m_I.setToIdentity();
    strain=(gradU+gradU.transpose())*0.5;// here the strain is small strain
}
void SmallStrainCahnHilliardMaterial::computeStressAndJacobian(const nlohmann::json &params,
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
//**************************************************************************
void SmallStrainCahnHilliardMaterial::computeFreeEnergyAndDerivatives(const nlohmann::json &parameters,
                                                 const VectorXd &args,
                                                 VectorXd       &F,
                                                 VectorXd       &dFdargs,
                                                 MatrixXd       &d2Fdargs2){
    double c,height;
    double ca,cb;
    c=args(1);
    height=JsonUtils::getValue(parameters,"Height");
    ca=JsonUtils::getValue(parameters,"Ca");
    cb=JsonUtils::getValue(parameters,"Cb");
    
    F(1)=height*(c-ca)*(c-ca)*(c-cb)*(c-cb);
    dFdargs(1)=height*2.0*(c-ca)*(c-cb)*(2*c-ca-cb);
    d2Fdargs2(1,1)=height*2.0*(ca*ca+4*ca*cb+cb*cb-6*ca*c-6*cb*c+6*c*c);

}