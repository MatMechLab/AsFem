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
//+++ Date   : 2022.08.22
//+++ Purpose: Implement the calculation of neo-hookean
//+++          hyperelastic material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/NeoHookeanMaterial.h"

void NeoHookeanMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||mate.getScalarMaterialsNum()){}
}

void NeoHookeanMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
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
    computeStressAndJacobian(inputparams,elmtinfo.m_dim,m_strain,m_pk2_stress,m_jacobian);

    m_stress=m_F*m_pk2_stress;// convert the 2nd PK stress to 1st PK stress
    m_devStress=m_stress.dev();
    m_devStrain=m_strain.dev();

    mate.ScalarMaterial("psi")=m_Psi;// for elastic free energy density

    mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devStress.doubledot(m_devStress));
    mate.ScalarMaterial("vonMises-strain")=sqrt(1.5*m_devStrain.doubledot(m_devStrain));
    mate.ScalarMaterial("hydrostatic-stress")=m_stress.trace()/3.0;

    mate.VectorMaterial("gradux")=elmtsoln.m_gpGradU[1];
    if(elmtinfo.m_dim>=2){
        mate.VectorMaterial("graduy")=elmtsoln.m_gpGradU[2];
        if(elmtinfo.m_dim==3) mate.VectorMaterial("graduz")=elmtsoln.m_gpGradU[3];
    }
    mate.Rank2Material("strain")=m_strain;
    mate.Rank2Material("stress")=m_stress;// the stress should be 1st PK stress, not 2nd PK stress !
    mate.Rank2Material("cauchy-stress")=m_F*m_pk2_stress*m_F.transpose()*(1.0/m_F.det());

    m_I.setToIdentity();
    mate.Rank4Material("jacobian")=m_I.ikXlj(m_pk2_stress)+m_jacobian.conjPushForward(m_F);// the final consistent jacobian

}

void NeoHookeanMaterial::computeStrain(const int &dim,const Rank2Tensor &gradU,Rank2Tensor &strain){
    if(dim){}
    m_I.setToIdentity();
    m_F=gradU+m_I;// deformation tensor F
    m_C=m_F.transpose()*m_F;// right Cauchy-Green strain C=F^tF
    m_Cinv=m_C.inverse();
    strain=(m_C-m_I)*0.5;// here the strain is E, the Lagrangian-Green strain
}
void NeoHookeanMaterial::computeStressAndJacobian(const nlohmann::json &params,
                                          const int &dim,
                                          const Rank2Tensor &strain,
                                          Rank2Tensor &stress,
                                          Rank4Tensor &jacobian){
    if(dim||strain.trace()){}
    double E=0.0,nu=0.0;
    double K=0.0,G=0.0;
    double lame=0.0;
    double J,I1;
    if(JsonUtils::hasValue(params,"E")&&
       JsonUtils::hasValue(params,"nu")){
        E=JsonUtils::getValue(params,"E");
        nu=JsonUtils::getValue(params,"nu");
        lame=E*nu/((1+nu)*(1-2*nu));
        G=0.5*E/(1.0+nu);
    }
    else if(JsonUtils::hasValue(params,"K")&&
            JsonUtils::hasValue(params,"G")){
        K=JsonUtils::getValue(params,"K");
        G=JsonUtils::getValue(params,"G");
        lame=K-G*2.0/3.0;
    }
    else if(JsonUtils::hasValue(params,"Lame")&&
            JsonUtils::hasValue(params,"mu")){
        lame=JsonUtils::getValue(params,"Lame");
        G=JsonUtils::getValue(params,"mu");
    }
    else{
        MessagePrinter::printErrorTxt("Invalid parameters, for neohookean material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }

    J=m_F.det();
    I1=m_C.trace();

    // here we use the polyconvex strain energy, it needs two parameters:
    // lame and G for the lame constants and shear moduli
    m_Psi=0.5*G*(I1-3.0)+(lame/4.0)*(J*J-1)-(0.5*lame+G)*log(J);
    // here the stress is 2nd PK stress, and the jacobian=dS/dE=2*dS/dC
    stress=m_Cinv*0.5*lame*(J*J-1.0)+(m_I-m_Cinv)*G;
    
    jacobian=m_Cinv.otimes(m_Cinv)*lame*J*J
            -m_Cinv.odot(m_Cinv)*lame*(J*J-1)
            +m_Cinv.odot(m_Cinv)*2.0*G;
    

}