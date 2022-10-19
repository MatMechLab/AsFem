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

#include "MateSystem/LinearElasticFractureMaterial.h"
#include "MathUtils/MathFuns.h"

LinearElasticFractureMaterial::LinearElasticFractureMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
LinearElasticFractureMaterial::~LinearElasticFractureMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
//******************************************************
void LinearElasticFractureMaterial::initMaterialProperties(const nlohmann::json &inputparams,
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
void LinearElasticFractureMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(mateold.getScalarMaterialsNum()){}

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
        MessagePrinter::printErrorTxt("LinearElasticFractureMaterial works only for 2d and 3d case, please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_I.setToIdentity();
    computeStrain(elmtinfo.m_dim,m_GradU,m_mechstrain);
    
    m_d=elmtsoln.m_gpU[1];

    m_args(1)=m_d;
    computeFreeEnergyAndDerivatives(inputparams,m_args,m_F,m_dFdargs,m_d2Fdargs2);
    
    computeStressAndJacobian(inputparams,elmtinfo.m_dim,m_mechstrain,m_stress,m_jacobian);

    m_devstress=m_stress.dev();

    mate.ScalarMaterial("F")=m_F(1);
    mate.ScalarMaterial("dFdD")=m_dFdargs(1);
    mate.ScalarMaterial("d2FdD2")=m_d2Fdargs2(1,1);

    mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devstress.doubledot(m_devstress));
    mate.ScalarMaterial("hydrostatic-stress")=m_stress.trace()/3.0;

    mate.Rank2Material("strain")=m_mechstrain;
    mate.Rank2Material("stress")=m_stress;
    mate.Rank2Material("dstressdD")=m_dstress_dD;

    mate.Rank4Material("jacobian")=m_jacobian;

    // for history variables
    if(m_psipos>mateold.ScalarMaterial("H")){
        mate.ScalarMaterial("H")=m_psipos;
        mate.Rank2Material("dHdstrain")=m_stress_pos;
    }
    else{
        mate.ScalarMaterial("H")=mateold.ScalarMaterial("H");
        mate.Rank2Material("dHdstrain").setToZeros();
    }

}
//**************************************************************************
void LinearElasticFractureMaterial::computeFreeEnergyAndDerivatives(const nlohmann::json &parameters,
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
void LinearElasticFractureMaterial::computeStrain(const int &dim,const Rank2Tensor &gradU,Rank2Tensor &strain){
    if(dim){}
    m_I.setToIdentity();
    strain=(gradU+gradU.transpose())*0.5;// here the strain is small strain
}
//**************************************************************************
void LinearElasticFractureMaterial::computeStressAndJacobian(const nlohmann::json &params,
                                          const int &dim,
                                          const Rank2Tensor &strain,
                                          Rank2Tensor &stress,
                                          Rank4Tensor &jacobian){
    if(dim){}

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
    else if(JsonUtils::hasValue(params,"Lame")&&
            JsonUtils::hasValue(params,"G")){
        lame=JsonUtils::getValue(params,"Lame");
        G=JsonUtils::getValue(params,"G");
    }
    else{
        MessagePrinter::printErrorTxt("Invalid parameters, for linear elastic fracture material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_projpos=strain.getPositiveProjectionTensor();
    m_i4sym.setToIdentity4Symmetric();
    m_projneg=m_i4sym-m_projpos;

    // strain splitting
    m_strain_pos=m_projpos.doubledot(strain);
    m_strain_neg=strain-m_strain_pos;

    double trEps,signpos,signneg;

    trEps=strain.trace();
    m_psipos=0.5*lame*MathFuns::bracketPos(trEps)*MathFuns::bracketPos(trEps)+G*(m_strain_pos*m_strain_pos).trace();
    m_psineg=0.5*lame*MathFuns::bracketNeg(trEps)*MathFuns::bracketNeg(trEps)+G*(m_strain_neg*m_strain_neg).trace();
    m_psi=g(m_d)*m_psipos+m_psineg;

    // for different stresses
    m_I.setToIdentity();// identity tensor
    m_stress_pos=m_I*lame*MathFuns::bracketPos(trEps)+m_strain_pos*2.0*G;
    m_stress_neg=m_I*lame*MathFuns::bracketNeg(trEps)+m_strain_neg*2.0*G;
    stress=m_stress_pos*(g(m_d)+m_stabilizer)+m_stress_neg;
    m_dstress_dD=m_stress_pos*dg(m_d);

    signpos=0.0;
    if(MathFuns::bracketPos(trEps)>0.0) signpos=1.0;

    signneg=0.0;
    if(MathFuns::bracketNeg(trEps)<0.0) signneg=1.0;

    jacobian=(m_I.otimes(m_I)*lame*signpos+m_projpos*2.0*G)*(g(m_d)+m_stabilizer)
             +m_I.otimes(m_I)*lame*signneg+m_projneg*2.0*G;

}