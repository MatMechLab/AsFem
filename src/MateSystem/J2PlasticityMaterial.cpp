//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.04.11
//+++ Purpose: Implement the J2 plasticity model
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/J2PlasticityMaterial.h"


void J2PlasticityMaterial::InitMaterialProperties(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln, Materials &Mate){
    // Here we do not consider any initial internal strains, stress
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}

    Mate.Rank2Materials("stress").SetToZeros();
    Mate.Rank2Materials("plastic_strain").SetToZeros();
    Mate.ScalarMaterials("effective_plastic_strain")=0;

}

//*********************************************************
void J2PlasticityMaterial::ComputeStrain(const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, RankTwoTensor &Strain){
    if(elmtinfo.nDim==1){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1]);
    }
    else if(elmtinfo.nDim==2){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2]);
    }
    else if(elmtinfo.nDim==3){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2],elmtsoln.gpGradU[3]);
    }
    Strain=(_GradU+_GradU.Transpose())*0.5;
}
//************************************************************
double J2PlasticityMaterial::ComputeYieldFunction(const vector<double> &InputParams, const double &trial_strss,
                                                  const double &effect_plastic_strain) {
    // parameters should be: E, nu, yield stress, hardening modulus
    const double YieldStress=InputParams[3-1];
    const double Hardening=InputParams[4-1];
    return trial_strss-sqrt(2.0/3.0)*(YieldStress+effect_plastic_strain*Hardening);
}
//****************************************************************
void J2PlasticityMaterial::ComputeAdmissibleStressState(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,const RankTwoTensor &Strain,const Materials &MateOld,Materials &Mate,RankTwoTensor &Stress,RankFourTensor &Jac){
    //************************************************
    // get rid of unused warnings
    //************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Strain(1,1)){}

    _E=InputParams[1-1];
    _nu=InputParams[2-1];
    _hardening_modulus=InputParams[4-1];

    // for bulk modulus and shear modulus
    _Lambda=_E/(3*(1-2*_nu)); // bulk modulus, not the first lame constant!!!
    _Mu=_E/(2*(1+_nu));       // shear modulus

    // for the variables in previous step
    _plastic_strain_old=MateOld.Rank2Materials("plastic_strain");
    _Effect_Plastic_Strain_Old=MateOld.ScalarMaterials("effective_plastic_strain");

    _I.SetToIdentity();
    _I4Sym.SetToIdentitySymmetric4();

    _devStrain=_Strain-_I*(Strain.Trace()/3.0);
    _STrial=(_devStrain-_plastic_strain_old)*2.0*_Mu;

    _F= ComputeYieldFunction(InputParams,_STrial.Norm(),_Effect_Plastic_Strain_Old);

    if(_F<=0.0){
        // for elastic case
        _N.SetToZeros();
        _DeltaGamma=0.0;
        _Jac.SetFromEandNu(_E,_nu);
    }
    else{
        // for plastic case
        _DeltaGamma=_F/(2*_Mu+2*_hardening_modulus/3.0);
        _N=_STrial/_STrial.Norm();
        _theta=1.0-2.0*_Mu*_DeltaGamma/_STrial.Norm();
        _thetabar=1.0/(1.0+_hardening_modulus/(3*_Mu))-(1-_theta);
        Jac=_I.OTimes(_I)*_Lambda
            +(_I4Sym-_I.OTimes(_I)*(1.0/3.0))*2*_Mu*_theta
            -_N.OTimes(_N)*2*_Mu*_thetabar;

    }
    // update all the variables
    Mate.ScalarMaterials("effective_plastic_strain")=_Effect_Plastic_Strain_Old+sqrt(2.0/3.0)*_DeltaGamma;
    Mate.Rank2Materials("plastic_strain")=_plastic_strain_old+_N*_DeltaGamma;
    Stress=_I*_Lambda*_Strain.Trace()+_STrial-_N*2*_Mu*_DeltaGamma;

}

void J2PlasticityMaterial::ComputeMaterialProperties(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,const Materials &MateOld,Materials &Mate) {

    //*********************************************************
    //*** get rid of unused warnings
    //*********************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}
 
    
    if(InputParams.size()<4){ 
        MessagePrinter::PrintErrorTxt("for J2-Plasticity material, four parameters are required, you need to give: E, nu, yield stress, and hardening modulus");
        MessagePrinter::AsFem_Exit();
    }


    ComputeStrain(elmtinfo,elmtsoln,_Strain);

    ComputeAdmissibleStressState(InputParams,elmtinfo,elmtsoln,_Strain,MateOld,Mate,_Stress,_Jac);

    Mate.Rank2Materials("stress")=_Stress;
    Mate.Rank2Materials("strain")=_Strain;
    Mate.Rank4Materials("jacobian")=_Jac;

    _devStress=_Stress-_I*(_Stress.Trace()/3.0);
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*_devStress.DoubleDot(_devStress));

}
