//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.04.11
//+++ Purpose: Implement the J2 plasticity model
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/J2PlasticityMaterial.h"

void J2PlasticityMaterial::InitMaterialProperties(const int &nDim, const Vector3d &gpCoord,
                                                  const vector<double> &InputParams, const vector<double> &gpU,
                                                  const vector<double> &gpUdot, const vector<Vector3d> &gpGradU,
                                                  const vector<Vector3d> &gpGradUdot, Materials &Mate) {
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUdot[0]||
       gpGradU[0](1)||gpGradUdot[0](1)||Mate.ScalarMaterials.size()){}

    Mate.Rank2Materials["plastic_strain"].SetToZeros();
    Mate.ScalarMaterials["effective_plastic_strain"]=0;
}
//********************************************************
void J2PlasticityMaterial::ComputeStrain(const int &nDim, const vector<Vector3d> &GradDisp, RankTwoTensor &Strain) {
    if(nDim==1){
        _GradU.SetFromGradU(GradDisp[1]);
    }
    else if(nDim==2){
        _GradU.SetFromGradU(GradDisp[1],GradDisp[2]);
    }
    else if(nDim==3){
        _GradU.SetFromGradU(GradDisp[1],GradDisp[2],GradDisp[3]);
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
void J2PlasticityMaterial::ComputeMaterialProperties(const double &t, const double &dt, const int &nDim,
                                                     const Vector3d &gpCoord, const vector<double> &InputParams,
                                                     const vector<double> &gpU, const vector<double> &gpUOld,
                                                     const vector<double> &gpUdot, const vector<double> &gpUdotOld,
                                                     const vector<Vector3d> &gpGradU,
                                                     const vector<Vector3d> &gpGradUOld,
                                                     const vector<Vector3d> &gpGradUdot,
                                                     const vector<Vector3d> &gpGradUdotOld, const Materials &MateOld,
                                                     Materials &Mate) {
    //**********************************************************************
    //*** get rid of unused warning
    //**********************************************************************
    if(t||dt||gpCoord(1)||gpU[0]||gpUOld[0]||
       gpUdot[0]||gpUdotOld[0]||gpGradU[0](1)||gpGradUOld[0](1)||
       gpGradUdot[0](1)||gpGradUdotOld[0](1)||MateOld.ScalarMaterials.size()){}

    if(InputParams.size()<4){
        MessagePrinter::PrintErrorTxt("for J2-Plasticity material, four parameters are required, you need to give: E, nu, yield stress, and hardening modulus");
        MessagePrinter::AsFem_Exit();
    }

    _E=InputParams[1-1];
    _nu=InputParams[2-1];
    _hardening_modulus=InputParams[4-1];

    // for bulk modulus and shear modulus
    _Lambda=_E/(3*(1-2*_nu)); // bulk modulus, not the first lame constant!!!
    _Mu=_E/(2*(1+_nu));       // shear modulus

    // for the variables in previous step
    _plastic_strain_old=MateOld.Rank2Materials.at("plastic_strain");
    _Effect_Plastic_Strain_Old=MateOld.ScalarMaterials.at("effective_plastic_strain");

    _I.SetToIdentity();
    _I4Sym.SetToIdentitySymmetric4();

    ComputeStrain(nDim,gpGradU,_Strain);

    _devStrain=_Strain-_I*(_Strain.Trace()/3.0);
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
        _Jac=_I.CrossDot(_I)*_Lambda
                +(_I4Sym-_I.CrossDot(_I)*(1.0/3.0))*2*_Mu*_theta
                -_N.CrossDot(_N)*2*_Mu*_thetabar;
    }
    // update all the variables
    Mate.ScalarMaterials["effective_plastic_strain"]=_Effect_Plastic_Strain_Old+sqrt(2.0/3.0)*_DeltaGamma;
    Mate.Rank2Materials["plastic_strain"]=_plastic_strain_old+_N*_DeltaGamma;
    Mate.Rank2Materials["stress"]=_I*_Lambda*_Strain.Trace()
            +_STrial-_N*2*_Mu*_DeltaGamma;
    _Stress=_I*_Lambda*_Strain.Trace()+_STrial-_N*2*_Mu*_DeltaGamma;
    Mate.Rank2Materials["strain"]=_Strain;
    Mate.Rank4Materials["jacobian"]=_Jac;

    _devStress=_Stress-_I*(_Stress.Trace()/3.0);
    Mate.ScalarMaterials["vonMises"]=sqrt(1.5*_devStress.DoubleDot(_devStress));

}