//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.11.25
//+++ Purpose: Calculate the material properties required by diffusion
//+++          coupled AC fracture element.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/DiffusionACFractureMaterial.h"
#include "MathUtils/MathFuns.h"

DiffusionACFractureMaterial::DiffusionACFractureMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
DiffusionACFractureMaterial::~DiffusionACFractureMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
//******************************************************
void DiffusionACFractureMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_Dt||elmtsoln.m_QpU[0]){}
    mate.ScalarMaterial("H")=0.0;
}
//********************************************************************
void DiffusionACFractureMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
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
        if(mate.BooleanMaterial("finite-strain")){
            MessagePrinter::printErrorTxt("Linear elastic fracture material works only in small strain case. Please check your input file and set 'finite-strain'=false");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        mate.BooleanMaterial("finite-strain")=false;// use small strain deformation as the default option
    }

    mate.ScalarMaterial("L")=JsonUtils::getValue(inputparams,"L");
    mate.ScalarMaterial("Gc")=JsonUtils::getValue(inputparams,"Gc");
    mate.ScalarMaterial("eps")=JsonUtils::getValue(inputparams,"eps");
    mate.ScalarMaterial("D")=JsonUtils::getValue(inputparams,"D");
    mate.ScalarMaterial("Omega")=JsonUtils::getValue(inputparams,"Omega");

    if(elmtinfo.m_Dim==2){
        m_GradU.setFromGradU2D(elmtsoln.m_QpGradU[3],elmtsoln.m_QpGradU[4]);// grad(ux), grad(uy)
    }
    else if(elmtinfo.m_Dim==3){
        m_GradU.setFromGradU3D(elmtsoln.m_QpGradU[3],elmtsoln.m_QpGradU[4],elmtsoln.m_QpGradU[5]);// grad(ux), grad(uy)
    }
    else{
        MessagePrinter::printErrorTxt("DiffusionACFractureMaterial works only for 2d and 3d case, please check your input file");
        MessagePrinter::exitAsFem();
    }

    m_I.setToIdentity();
    computeStrain(elmtinfo.m_Dim,m_GradU,m_strain);
    
    m_d=elmtsoln.m_QpU[2];

    m_args(1)=m_d;
    computeFreeEnergyAndDerivatives(inputparams,m_args,m_F,m_dFdargs,m_d2Fdargs2);

    m_args(1)=elmtsoln.m_QpU[1];// stores the concentration
    
    computeStressAndJacobian(inputparams,elmtinfo.m_Dim,m_strain,m_stress,m_jacobian);

    if(elmtinfo.m_Dim==2){
        m_strain(3,3)=m_eps_zz;
    }

    m_devstress=m_stress.dev();

    mate.ScalarMaterial("F")=m_F(1);
    mate.ScalarMaterial("dFdD")=m_dFdargs(1);
    mate.ScalarMaterial("d2FdD2")=m_d2Fdargs2(1,1);

    mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devstress.doubledot(m_devstress));
    mate.ScalarMaterial("hydrostatic-stress")=m_stress.trace()/3.0;

    mate.Rank2Material("strain")=m_strain;
    mate.Rank2Material("stress")=m_stress;
    mate.Rank2Material("dstressdD")=m_dstress_dD;
    mate.Rank2Material("dstressdc")=m_dstress_dc;

    mate.ScalarMaterial("SigmaH")=m_dstress_dc.trace()/3.0;// this is the factor of hydrostatic stress, not itself !!!

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
void DiffusionACFractureMaterial::computeFreeEnergyAndDerivatives(const nlohmann::json &parameters,
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
void DiffusionACFractureMaterial::computeStrain(const int &dim,const Rank2Tensor &gradU,Rank2Tensor &strain){
    if(dim){}
    m_I.setToIdentity();
    strain=(gradU+gradU.transpose())*0.5;// here the strain is small strain
}
//**************************************************************************
void DiffusionACFractureMaterial::computeStressAndJacobian(const nlohmann::json &params,
                                          const int &dim,
                                          const Rank2Tensor &strain,
                                          Rank2Tensor &stress,
                                          Rank4Tensor &jacobian){

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
        E=9.0*K*G/(3.0*K+G);
        nu=(3.0*K-2.0*G)/(2.0*(3.0*K+G));
    }
    else if(JsonUtils::hasValue(params,"Lame")&&
            JsonUtils::hasValue(params,"mu")){
        lame=JsonUtils::getValue(params,"Lame");
        G=JsonUtils::getValue(params,"mu");
        K=lame+2.0*G/3.0;
        E=9.0*K*G/(3.0*K+G);
        nu=0.5*(3.0*K-2.0*G)/(3.0*K+G);
    }
    else if(JsonUtils::hasValue(params,"Lame")&&
            JsonUtils::hasValue(params,"G")){
        lame=JsonUtils::getValue(params,"Lame");
        G=JsonUtils::getValue(params,"G");
        K=lame+2.0*G/3.0;
        E=G*(3.0*lame+2.0*G)/(lame+G);
        nu=0.5*lame/(lame+G);
    }
    else{
        MessagePrinter::printErrorTxt("Invalid parameters, for diffusion allencahn fracture material, you should give either E,nu or K,G or Lame,G. Please check your input file");
        MessagePrinter::exitAsFem();
    }

    double trEps,dtrEpsdc,signpos,signneg;

    m_sig_zz=0.0;m_eps_zz=0.0;

    m_I.setToIdentity();
    // for the concentration induced eigen strain
    m_Omega=JsonUtils::getValue(params,"Omega");
    m_Cref=JsonUtils::getValue(params,"Cref");

    m_mech_strain=strain-m_I*(m_args(1)-m_Cref)*m_Omega/3.0;
    m_dmechstrain_dc=m_I*-1.0*m_Omega/3.0;

    if(dim!=2){
        m_projpos=m_mech_strain.getPositiveProjectionTensor();
        m_i4sym.setToIdentity4Symmetric();
        m_projneg=m_i4sym-m_projpos;

        // strain splitting
        m_strain_pos=m_projpos.doubledot(m_mech_strain);
        m_strain_neg=m_mech_strain-m_strain_pos;

        m_dstrainpos_dc=m_projpos.doubledot(m_dmechstrain_dc);
        m_dstrainneg_dc=m_dmechstrain_dc-m_dstrainpos_dc;

        trEps=m_mech_strain.trace();
        dtrEpsdc=m_dmechstrain_dc.trace();
        m_psipos=0.5*lame*MathFuns::bracketPos(trEps)*MathFuns::bracketPos(trEps)+G*(m_strain_pos*m_strain_pos).trace();
        m_psineg=0.5*lame*MathFuns::bracketNeg(trEps)*MathFuns::bracketNeg(trEps)+G*(m_strain_neg*m_strain_neg).trace();
        m_psi=g(m_d)*m_psipos+m_psineg;

        // for different stresses
        m_I.setToIdentity();// identity tensor
        m_stress_pos=m_I*lame*MathFuns::bracketPos(trEps)+m_strain_pos*2.0*G;
        m_stress_neg=m_I*lame*MathFuns::bracketNeg(trEps)+m_strain_neg*2.0*G;

        m_dstresspos_dc=m_I*lame*MathFuns::bracketPos(dtrEpsdc)+m_dstrainpos_dc*2.0*G;
        m_dstressneg_dc=m_I*lame*MathFuns::bracketNeg(dtrEpsdc)+m_dstrainneg_dc*2.0*G;

        stress=m_stress_pos*(g(m_d)+m_stabilizer)+m_stress_neg;
        m_dstress_dD=m_stress_pos*dg(m_d);
        m_dstress_dc=m_dstresspos_dc*(g(m_d)+m_stabilizer)+m_dstressneg_dc;

        signpos=0.0;
        if(MathFuns::bracketPos(trEps)>0.0) signpos=1.0;

        signneg=0.0;
        if(MathFuns::bracketNeg(trEps)<0.0) signneg=1.0;

        jacobian=(m_I.otimes(m_I)*lame*signpos+m_projpos*2.0*G)*(g(m_d)+m_stabilizer)
                 +m_I.otimes(m_I)*lame*signneg+m_projneg*2.0*G;
    }
    else{
        // for 2d case
        if(JsonUtils::hasValue(params,"plane-strain")){
            if(JsonUtils::getBoolean(params,"plane-strain")){
                // for plane strain case
                m_sig_zz=(E/(1.0+nu))*(nu/(1.0-2.0*nu))*(strain(1,1)+strain(2,2));
                m_eps_zz=0.0;
                
                E=E/(1.0-nu*nu);
                nu=nu/(1.0-nu);

                lame=E*nu/((1+nu)*(1-2*nu));
                G=0.5*E/(1.0+nu);

                m_projpos=m_mech_strain.getPositiveProjectionTensor();
                m_i4sym.setToIdentity4Symmetric();
                m_projneg=m_i4sym-m_projpos;

                // strain splitting
                m_strain_pos=m_projpos.doubledot(m_mech_strain);
                m_strain_neg=m_mech_strain-m_strain_pos;
                m_dstrainpos_dc=m_projpos.doubledot(m_dmechstrain_dc);
                m_dstrainneg_dc=m_dmechstrain_dc-m_dstrainpos_dc;

                trEps=m_mech_strain.trace();
                dtrEpsdc=m_dmechstrain_dc.trace();
                m_psipos=0.5*lame*MathFuns::bracketPos(trEps)*MathFuns::bracketPos(trEps)+G*(m_strain_pos*m_strain_pos).trace();
                m_psineg=0.5*lame*MathFuns::bracketNeg(trEps)*MathFuns::bracketNeg(trEps)+G*(m_strain_neg*m_strain_neg).trace();
                m_psi=g(m_d)*m_psipos+m_psineg;

                // for different stresses
                m_I.setToIdentity();// identity tensor
                m_stress_pos=m_I*lame*MathFuns::bracketPos(trEps)+m_strain_pos*2.0*G;
                m_stress_neg=m_I*lame*MathFuns::bracketNeg(trEps)+m_strain_neg*2.0*G;
                m_dstresspos_dc=m_I*lame*MathFuns::bracketPos(dtrEpsdc)+m_dstrainpos_dc*2.0*G;
                m_dstressneg_dc=m_I*lame*MathFuns::bracketNeg(dtrEpsdc)+m_dstrainneg_dc*2.0*G;
                stress=m_stress_pos*(g(m_d)+m_stabilizer)+m_stress_neg;
                stress(3,3)=m_sig_zz;
                m_dstress_dD=m_stress_pos*dg(m_d);
                m_dstress_dc=m_dstresspos_dc*(g(m_d)+m_stabilizer)+m_dstressneg_dc;

                signpos=0.0;
                if(MathFuns::bracketPos(trEps)>0.0) signpos=1.0;

                signneg=0.0;
                if(MathFuns::bracketNeg(trEps)<0.0) signneg=1.0;
                
                jacobian=(m_I.otimes(m_I)*lame*signpos+m_projpos*2.0*G)*(g(m_d)+m_stabilizer)
                         +m_I.otimes(m_I)*lame*signneg+m_projneg*2.0*G;
            }
            else{
                // for plane-stress case
                m_projpos=m_mech_strain.getPositiveProjectionTensor();
                m_i4sym.setToIdentity4Symmetric();
                m_projneg=m_i4sym-m_projpos;

                // strain splitting
                m_strain_pos=m_projpos.doubledot(m_mech_strain);
                m_strain_neg=m_mech_strain-m_strain_pos;
                m_dstrainpos_dc=m_projpos.doubledot(m_dmechstrain_dc);
                m_dstrainneg_dc=m_dmechstrain_dc-m_dstrainpos_dc;

                trEps=m_mech_strain.trace();
                dtrEpsdc=m_dmechstrain_dc.trace();
                m_psipos=0.5*lame*MathFuns::bracketPos(trEps)*MathFuns::bracketPos(trEps)+G*(m_strain_pos*m_strain_pos).trace();
                m_psineg=0.5*lame*MathFuns::bracketNeg(trEps)*MathFuns::bracketNeg(trEps)+G*(m_strain_neg*m_strain_neg).trace();
                m_psi=g(m_d)*m_psipos+m_psineg;

                // for different stresses
                m_I.setToIdentity();// identity tensor
                m_stress_pos=m_I*lame*MathFuns::bracketPos(trEps)+m_strain_pos*2.0*G;
                m_stress_neg=m_I*lame*MathFuns::bracketNeg(trEps)+m_strain_neg*2.0*G;
                m_dstresspos_dc=m_I*lame*MathFuns::bracketPos(dtrEpsdc)+m_dstrainpos_dc*2.0*G;
                m_dstressneg_dc=m_I*lame*MathFuns::bracketNeg(dtrEpsdc)+m_dstrainneg_dc*2.0*G;
                stress=m_stress_pos*(g(m_d)+m_stabilizer)+m_stress_neg;
                m_dstress_dD=m_stress_pos*dg(m_d);
                m_dstress_dc=m_dstresspos_dc*(g(m_d)+m_stabilizer)+m_dstressneg_dc;

                signpos=0.0;
                if(MathFuns::bracketPos(trEps)>0.0) signpos=1.0;

                signneg=0.0;
                if(MathFuns::bracketNeg(trEps)<0.0) signneg=1.0;
                
                jacobian=(m_I.otimes(m_I)*lame*signpos+m_projpos*2.0*G)*(g(m_d)+m_stabilizer)
                         +m_I.otimes(m_I)*lame*signneg+m_projneg*2.0*G;

                m_eps_zz=(-nu/E)*(stress(1,1)+stress(2,2));
            }
        }
        else{
            // default option is plane-stress
            m_projpos=m_mech_strain.getPositiveProjectionTensor();
            m_i4sym.setToIdentity4Symmetric();
            m_projneg=m_i4sym-m_projpos;

            // strain splitting
            m_strain_pos=m_projpos.doubledot(m_mech_strain);
            m_strain_neg=m_mech_strain-m_strain_pos;
            m_dstrainpos_dc=m_projpos.doubledot(m_dmechstrain_dc);
            m_dstrainneg_dc=m_dmechstrain_dc-m_dstrainpos_dc;

            trEps=m_mech_strain.trace();
            dtrEpsdc=m_dmechstrain_dc.trace();
            m_psipos=0.5*lame*MathFuns::bracketPos(trEps)*MathFuns::bracketPos(trEps)+G*(m_strain_pos*m_strain_pos).trace();
            m_psineg=0.5*lame*MathFuns::bracketNeg(trEps)*MathFuns::bracketNeg(trEps)+G*(m_strain_neg*m_strain_neg).trace();
            m_psi=g(m_d)*m_psipos+m_psineg;

            // for different stresses
            m_I.setToIdentity();// identity tensor
            m_stress_pos=m_I*lame*MathFuns::bracketPos(trEps)+m_strain_pos*2.0*G;
            m_stress_neg=m_I*lame*MathFuns::bracketNeg(trEps)+m_strain_neg*2.0*G;
            m_dstresspos_dc=m_I*lame*MathFuns::bracketPos(dtrEpsdc)+m_dstrainpos_dc*2.0*G;
            m_dstressneg_dc=m_I*lame*MathFuns::bracketNeg(dtrEpsdc)+m_dstrainneg_dc*2.0*G;
            stress=m_stress_pos*(g(m_d)+m_stabilizer)+m_stress_neg;
            m_dstress_dD=m_stress_pos*dg(m_d);
            m_dstress_dc=m_dstresspos_dc*(g(m_d)+m_stabilizer)+m_dstressneg_dc;

            signpos=0.0;
            if(MathFuns::bracketPos(trEps)>0.0) signpos=1.0;

            signneg=0.0;
            if(MathFuns::bracketNeg(trEps)<0.0) signneg=1.0;
                
            jacobian=(m_I.otimes(m_I)*lame*signpos+m_projpos*2.0*G)*(g(m_d)+m_stabilizer)
                     +m_I.otimes(m_I)*lame*signneg+m_projneg*2.0*G;

            m_eps_zz=(-nu/E)*(stress(1,1)+stress(2,2));
        }
    }
}