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
//+++ Purpose: Implement 1d plasticity model with linear hardening
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/Plastic1DMaterial.h"

void Plastic1DMaterial::InitMaterialProperties(const int &nDim, const Vector3d &gpCoord,
                                               const vector<double> &InputParams, const vector<double> &gpU,
                                               const vector<double> &gpUdot, const vector<Vector3d> &gpGradU,
                                               const vector<Vector3d> &gpGradUdot, Materials &Mate) {
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUdot[0]||
       gpGradU[0](1)||gpGradUdot[0](1)||Mate.ScalarMaterials.size()){}

    Mate.ScalarMaterials["effective_plastic_strain"]=0;
    Mate.Rank2Materials["plastic_strain"].SetToZeros();
}
//****************************************************************************
void Plastic1DMaterial::ComputeStrain(const int &nDim, const vector<Vector3d> &GradDisp, RankTwoTensor &Strain) {
    if(nDim==1){
        _GradU.SetFromGradU(GradDisp[1]);
    }
    else {
        MessagePrinter::PrintErrorTxt("Plastic1DMaterial only works for 1d case");
        MessagePrinter::AsFem_Exit();
    }
    Strain=(_GradU+_GradU.Transpose())*0.5;
}
//****************************************************************************
double Plastic1DMaterial::ComputeYieldFunction(const vector<double> &InputParams, const double &trial_strss,
                                               const double &effect_plastic_strain) {
    // parameters should be: E, nu, yield stress, hardening modulus
    const double YieldStress=InputParams[2-1];
    const double Hardening=InputParams[3-1];
    return trial_strss-(YieldStress+effect_plastic_strain*Hardening);
}
//****************************************************************************
void Plastic1DMaterial::ComputeMaterialProperties(const double &t, const double &dt, const int &nDim,
                                                  const Vector3d &gpCoord, const vector<double> &InputParams,
                                                  const vector<double> &gpU, const vector<double> &gpUOld,
                                                  const vector<double> &gpUdot, const vector<double> &gpUdotOld,
                                                  const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradUOld,
                                                  const vector<Vector3d> &gpGradUdot,
                                                  const vector<Vector3d> &gpGradUdotOld, const Materials &MateOld,
                                                  Materials &Mate) {
    //**********************************************************************
    //*** get rid of unused warning
    //**********************************************************************
    if(t||dt||gpCoord(1)||gpU[0]||gpUOld[0]||
       gpUdot[0]||gpUdotOld[0]||gpGradU[0](1)||gpGradUOld[0](1)||
       gpGradUdot[0](1)||gpGradUdotOld[0](1)||MateOld.ScalarMaterials.size()){}

    if(InputParams.size()<3){
        MessagePrinter::PrintErrorTxt("for 1d plastic material, three parameters are required, you need to give: E, yield stress, and hardening modulus");
        MessagePrinter::AsFem_Exit();
    }

    _E=InputParams[1-1];
    _hardening_modulus=InputParams[3-1];

    ComputeStrain(nDim,gpGradU,_Strain);// this is the total strain

    _PlasticStrainOld=MateOld.ScalarMaterials.at("effective_plastic_strain");
    _TrialStrain=_Strain(1,1)-MateOld.Rank2Materials.at("plastic_strain")(1,1);
    _TrialStress=_E*_TrialStrain;

    _F= ComputeYieldFunction(InputParams,_TrialStress,_PlasticStrainOld);
    if(_F<=0.0){
        // for elastic loading/unloading
        _Jac.SetToZeros();
        _Jac(1,1,1,1)=_E;
        _DeltaGamma=0.0;
    }
    else{
        // for plastic deformation
        _DeltaGamma=_F/(_E+_hardening_modulus);
        _Jac.SetToZeros();
        _Jac(1,1,1,1)=_E*_hardening_modulus/(_E+_hardening_modulus);
    }

    _Stress.SetToZeros();
    _Stress(1,1)=_TrialStress-_E*_DeltaGamma*Sign(_TrialStress);

    Mate.Rank2Materials["plastic_strain"].SetToZeros();
    Mate.Rank2Materials["plastic_strain"](1,1)=MateOld.Rank2Materials.at("plastic_strain")(1,1)+_DeltaGamma*Sign(_TrialStress);

    Mate.ScalarMaterials["effective_plastic_strain"]=MateOld.ScalarMaterials.at("effective_plastic_strain")+_DeltaGamma;

    _I.SetToIdentity();
    _devStress=_Stress-_I*(_Stress.Trace()/3.0);
    Mate.ScalarMaterials["vonMises"]=sqrt(1.5*_devStress.DoubleDot(_devStress));

    Mate.Rank2Materials["stress"]=_Stress;
    Mate.Rank2Materials["strain"]=_Strain;
    Mate.Rank4Materials["jacobian"]=_Jac;

}