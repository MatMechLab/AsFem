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

#include "MateSystem/NeoHookeanPFFractureMaterial.h"
#include "MathUtils/MathFuns.h"

NeoHookeanPFFractureMaterial::NeoHookeanPFFractureMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
NeoHookeanPFFractureMaterial::~NeoHookeanPFFractureMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
//******************************************************
void NeoHookeanPFFractureMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]){}
    mate.ScalarMaterial("H")=0.0;
}
//********************************************************************
void NeoHookeanPFFractureMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
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
        if(!mate.BooleanMaterial("finite-strain")){
            MessagePrinter::printErrorTxt("For neohookean fracture material, you must enable finite-strain option. Please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        mate.BooleanMaterial("finite-strain")=true;// use finite strain deformation as the default option
    }

    mate.ScalarMaterial("L")=JsonUtils::getValue(inputparams,"L");
    mate.ScalarMaterial("Gc")=JsonUtils::getValue(inputparams,"Gc");
    mate.ScalarMaterial("eps")=JsonUtils::getValue(inputparams,"eps");

    if(elmtinfo.m_dim==2){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[2],elmtsoln.m_gpGradU[3]);// grad(ux), grad(uy)
    }
    else if(elmtinfo.m_dim==3){
        m_GradU.setFromGradU(elmtsoln.m_gpGradU[2],elmtsoln.m_gpGradU[3],elmtsoln.m_gpGradU[4]);// grad(ux), grad(uy)
    }
    else{
        MessagePrinter::printErrorTxt("NeoHookeanPFFractureMaterial works only for 2d and 3d case, please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_I.setToIdentity();
    computeStrain(elmtinfo.m_dim,m_GradU,m_Estrain);
    
    m_d=elmtsoln.m_gpU[1];

    m_args(1)=m_d;
    computeFreeEnergyAndDerivatives(inputparams,m_args,m_F,m_dFdargs,m_d2Fdargs2);
    
    computeStressAndJacobian(inputparams,elmtinfo.m_dim,m_Estrain,m_PK1stress,m_jacobian);

    m_devstress=m_PK1stress.dev();

    mate.ScalarMaterial("F")=m_F(1);
    mate.ScalarMaterial("dFdD")=m_dFdargs(1);
    mate.ScalarMaterial("d2FdD2")=m_d2Fdargs2(1,1);

    mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devstress.doubledot(m_devstress));
    mate.ScalarMaterial("hydrostatic-stress")=m_PK1stress.trace()/3.0;

    mate.Rank2Material("strain")=m_Estrain;
    mate.Rank2Material("stress")=m_PK1stress;
    mate.Rank2Material("dstressdD")=m_dPK1stress_dD;

    mate.Rank4Material("jacobian")=m_I.ikXlj(m_PK2stress)+m_jacobian.conjPushForward(m_Fe);

    // for history variables
    if(m_psipos>mateold.ScalarMaterial("H")){
        mate.ScalarMaterial("H")=m_psipos;
        mate.Rank2Material("dHdstrain")=m_Fe*m_PK2stress_pos;
    }
    else{
        mate.ScalarMaterial("H")=mateold.ScalarMaterial("H");
        mate.Rank2Material("dHdstrain").setToZeros();
    }

}
//**************************************************************************
void NeoHookeanPFFractureMaterial::computeFreeEnergyAndDerivatives(const nlohmann::json &parameters,
                                                 const VectorXd &args,
                                                 VectorXd       &F,
                                                 VectorXd       &dFdargs,
                                                 MatrixXd       &d2Fdargs2){
    double d,Gc,eps;
    d=args(1);
    Gc=JsonUtils::getValue(parameters,"Gc");
    eps=JsonUtils::getValue(parameters,"eps");
    F(1)=0.5*d*d*Gc/eps;
    dFdargs(1)=d*Gc/eps;
    d2Fdargs2(1,1)=Gc/eps;
}
//**************************************************************************
void NeoHookeanPFFractureMaterial::computeStrain(const int &dim,const Rank2Tensor &gradU,Rank2Tensor &strain){
    if(dim){}
    m_I.setToIdentity();
    m_Fe=gradU+m_I;
    m_Ce=m_Fe.transpose()*m_Fe;// Ce=Fe^t*Fe
    strain=(m_Ce-m_I)*0.5;// here the strain is E, the Lagrangian-Green strain
}
//**************************************************************************
void NeoHookeanPFFractureMaterial::computeStressAndJacobian(const nlohmann::json &params,
                                          const int &dim,
                                          const Rank2Tensor &strain,
                                          Rank2Tensor &stress,
                                          Rank4Tensor &jacobian){
    if(dim||strain(1,1)){}

    m_stabilizer=JsonUtils::getValue(params,"stabilizer");
    double E=0.0,nu=0.0;
    double K=0.0,G=0.0;
    double lame=0.0;

    jacobian.setToZeros();

    if(JsonUtils::hasValue(params,"E")&&
       JsonUtils::hasValue(params,"nu")){
        E=JsonUtils::getValue(params,"E");
        nu=JsonUtils::getValue(params,"nu");
        lame=E*nu/((1+nu)*(1-2*nu));
        K=E/(3.0*(1.0-2.0*nu));
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
        K=lame+2.0*G/3.0;
    }
    else if(JsonUtils::hasValue(params,"Lame")&&
            JsonUtils::hasValue(params,"G")){
        lame=JsonUtils::getValue(params,"Lame");
        G=JsonUtils::getValue(params,"G");
        K=lame+2.0*G/3.0;
    }
    else{
        MessagePrinter::printErrorTxt("Invalid parameters, for neohookean pf fracture material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_Je=m_Fe.det();
    m_Je23=std::pow(m_Je,-2.0/3.0);
    m_I1=m_Ce.trace();
    m_I1bar=m_I1*m_Je23;

    m_I.setToIdentity();// identity tensor
    m_CeInv=m_Ce.inverse();// inverse of Ce

    if(m_Je>1.0){
        // for tensile loading
        m_psipos=0.5*K*(0.5*(m_Je*m_Je-1.0)-std::log(m_Je))
                +0.5*G*(m_I1bar-3.0);
        m_psineg=0.0;
        
        m_PK2stress_pos=m_CeInv*0.5*K*(m_Je*m_Je-1.0)
                       -m_CeInv*(G/3.0)*m_Je23*m_I1
                       +m_I*G*m_Je23;
        m_PK2stress_neg.setToZeros();

        m_jacobian_pos=m_CeInv.otimes(m_CeInv)*K*m_Je*m_Je
                      -m_CeInv.odot(m_CeInv)*K*(m_Je*m_Je-1.0)
                      +(
                        m_CeInv.otimes(m_CeInv)*(m_I1/3.0)
                       -m_CeInv.otimes(m_I)
                       +m_CeInv.odot(m_CeInv)*m_I1
                       -m_I.otimes(m_CeInv)
                      )*(2.0*G/3.0)*m_Je23;
        m_jacobian_neg.setToZeros();
    }
    else{
        // for compressive case
        m_psipos=0.5*G*(m_I1bar-3.0);
        m_psineg=0.5*K*(0.5*(m_Je*m_Je-1.0)-std::log(m_Je));

        m_PK2stress_pos=m_I*G*m_Je23-m_CeInv*(G/3.0)*m_Je23*m_I1;
        m_PK2stress_neg=m_CeInv*K*0.5*(m_Je*m_Je-1.0);

        m_jacobian_pos=(
                        m_CeInv.otimes(m_CeInv)*(m_I1/3.0)
                       -m_CeInv.otimes(m_I)
                       +m_CeInv.odot(m_CeInv)*m_I1
                       -m_I.otimes(m_CeInv)
                      )*(2.0*G/3.0)*m_Je23;
        m_jacobian_neg=m_CeInv.otimes(m_CeInv)*K*m_Je*m_Je
                      -m_CeInv.odot(m_CeInv)*K*(m_Je*m_Je-1.0);
    }

    m_psi=m_psipos*(g(m_d)+m_stabilizer)+m_psineg;

    // for different stresses and their derivatives
    m_PK2stress=m_PK2stress_pos*(g(m_d)+m_stabilizer)+m_PK2stress_neg;
    m_dPK2stress_dD=m_PK2stress_pos*dg(m_d);

    stress=m_Fe*m_PK2stress;// the 1st PK stress
    m_dPK1stress_dD=m_Fe*m_dPK2stress_dD;// the 1st PK stress's derivative wrt d

    // for the jacobian in reference configuration
    jacobian=m_jacobian_pos*(g(m_d)+m_stabilizer)+m_jacobian_neg;

    // TODO: add plane-stress modification !

}