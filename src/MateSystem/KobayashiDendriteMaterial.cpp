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
//+++ Date   : 2022.11.12
//+++ Purpose: Calculate the free energy and its derivatives based
//+++          Kobayashi's dendrite model
//+++ Ref    : Modeling and numerical simulations of dendritic crystal growth
//+++ DOI    : https://doi.org/10.1016/0167-2789(93)90120-P
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "MateSystem/KobayashiDendriteMaterial.h"

KobayashiDendriteMaterial::KobayashiDendriteMaterial(){
    m_args.resize(11);
    m_F.resize(11);
    m_dFdargs.resize(11);
    m_d2Fdargs2.resize(11,11);
}
KobayashiDendriteMaterial::~KobayashiDendriteMaterial(){
    m_args.clean();
    m_F.clean();
    m_dFdargs.clean();
    m_d2Fdargs2.clean();
}

void KobayashiDendriteMaterial::initMaterialProperties(const nlohmann::json &inputparams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        MaterialsContainer &mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(inputparams.size()||elmtinfo.m_dt||elmtsoln.m_gpU[0]||mate.getScalarMaterialsNum()){}

}

//********************************************************************
void KobayashiDendriteMaterial::computeMaterialProperties(const nlohmann::json &inputparams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const MaterialsContainer &mateold,
                                           MaterialsContainer &mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(mateold.getScalarMaterialsNum()){}

    if(elmtinfo.m_dim!=2){
        MessagePrinter::printErrorTxt("KobayashiDendrite material only works for 2d case");
        MessagePrinter::exitAsFem();
    }

    if(!JsonUtils::hasOnlyGivenValues(inputparams,vector<string>{"L","k0","delta","N","Latent-heat"})){
        MessagePrinter::printErrorTxt("KobayashiDendrite material requires only: L, k0, delta, N, and Latent-heat "
                                      "please check your input file");
        MessagePrinter::exitAsFem();
    }

    mate.ScalarMaterial("L")=JsonUtils::getValue(inputparams,"L");
    mate.ScalarMaterial("Latent-heat")=JsonUtils::getValue(inputparams,"Latent-heat");
    
    m_args(1)=elmtsoln.m_gpU[1];// for order parameter
    m_args(2)=elmtsoln.m_gpU[2];// for temperature
    computeFreeEnergyAndDerivatives(inputparams,m_args,m_F,m_dFdargs,m_d2Fdargs2);

    mate.ScalarMaterial("F")=m_F(1);
    // for the derivatives w.r.t. eta
    mate.ScalarMaterial("dFdeta")=m_dFdargs(1);
    mate.ScalarMaterial("d2Fdeta2")=m_d2Fdargs2(1,1);
    mate.ScalarMaterial("d2FdetadT")=m_d2Fdargs2(1,2);

    // for the anisotropic gradient and its derivatives
    double K0,K,dK,ddK,N,delta;
    const double tol=1.0e-9;
    double threshold;
    double norm,normsq,n;
    Vector3d GradEta,dKdGradEta,ddKdGradEta;

    K0=JsonUtils::getValue(inputparams,"k0");
    N=JsonUtils::getValue(inputparams,"N");
    delta=JsonUtils::getValue(inputparams,"delta");

    GradEta=elmtsoln.m_gpGradU[1];
    threshold=1.0-tol;
    normsq=GradEta(1)*GradEta(1)+GradEta(2)*GradEta(2);
    norm=sqrt(normsq);

    n=0.0;
    if(normsq>tol){
        n=GradEta(1)/norm;
    }

    if(n> threshold) n= threshold;
    if(n<-threshold) n=-threshold;

    double theta,dthetadn;
    Vector3d dndgradeta;

    theta   = std::acos(n)*MathFuns::sign(GradEta(2));
    dthetadn=-MathFuns::sign(GradEta(2))/std::sqrt(1.0-n*n);

    dndgradeta=0.0;
    if(normsq>tol){
        dndgradeta(1)= GradEta(2)*GradEta(2)/(norm*normsq);
        dndgradeta(2)=-GradEta(1)*GradEta(2)/(norm*normsq);
        dndgradeta(3)=0.0;
    }

    K  = K0*(1.0+delta*cos(N*theta));
    dK =-K0*delta*N*sin(N*theta);
    ddK=-K0*delta*N*N*cos(N*theta);

    dKdGradEta=dK*dthetadn*dndgradeta;
    ddKdGradEta=ddK*dthetadn*dndgradeta;

    // store the variables into different materials
    mate.ScalarMaterial("K")=K;
    mate.ScalarMaterial("dK")=dK;
    mate.VectorMaterial("dKdGradEta")=dKdGradEta;
    mate.VectorMaterial("ddKdGradEta")=ddKdGradEta;
    
    mate.VectorMaterial("gradeta")=elmtsoln.m_gpGradU[1];// the gradient of eta
    mate.VectorMaterial("gradT")  =elmtsoln.m_gpGradU[2];// the gradient of T

}
//**************************************************************
void KobayashiDendriteMaterial::computeFreeEnergyAndDerivatives(const nlohmann::json &parameters,
                                                                  const VectorXd &args,
                                                                  VectorXd       &F,
                                                                  VectorXd       &dFdargs,
                                                                  MatrixXd       &d2Fdargs2){
    if(parameters.size()){}
    double eta=args(1);
    double T  =args(2);
    const double PI=3.141592653589793238462;
    double m,dmdT;
    // Now we start to calculate the double well free energy density(T-dependent!)
    m=0.9*atan(10.0*T)/PI;
    dmdT=0.09/(PI*T*T+0.01*PI);
    
    F(1)=0.25*pow(eta,4)-(0.5-m/3.0)*pow(eta,3)+(0.25-0.5*m)*eta*eta;//F
    dFdargs(1)=0.5*(eta-1.0)*eta*(2*eta+2*m-1.0);//dF/deta
    d2Fdargs2(1,1)=m*(2*eta-1.0)+3*eta*eta-3*eta+0.5;//d^2F/deta^2
    d2Fdargs2(1,2)=0.5*(eta-1.0)*eta*2*dmdT;// d^2F/(detadT)

}