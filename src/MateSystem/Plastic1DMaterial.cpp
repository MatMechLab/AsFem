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


void Plastic1DMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    // Here we do not consider any initial internal strains, stress
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}


    Mate.Rank2Materials("plastic_strain").SetToZeros();
    Mate.ScalarMaterials("effective_plastic_strain")=0;

}
//*********************************************************
void Plastic1DMaterial::ComputeStrain(const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, RankTwoTensor &Strain){
    if(elmtinfo.nDim==1){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1]);
    }
    else{
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
void Plastic1DMaterial::ComputeAdmissibleStressState(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,const RankTwoTensor &Strain,const Materials &MateOld,Materials &Mate,RankTwoTensor &Stress,RankFourTensor &Jac){
    //************************************************
    // get rid of unused warnings
    //************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Strain(1,1)){}

    _E=InputParams[1-1];
    _hardening_modulus=InputParams[4-1];
    
    _PlasticStrainOld=MateOld.ScalarMaterials("effective_plastic_strain");
    _TrialStrain=_Strain(1,1)-MateOld.Rank2Materials("plastic_strain")(1,1);
    _TrialStress=_E*_TrialStrain;

    _F= ComputeYieldFunction(InputParams,_TrialStress,_PlasticStrainOld);
    if(_F<=0.0){
        // for elastic loading/unloading
        Jac.SetToZeros();
        Jac(1,1,1,1)=_E;
        _DeltaGamma=0.0;
    }
    else{
        // for plastic deformation
        _DeltaGamma=_F/(_E+_hardening_modulus);
        Jac.SetToZeros();
        Jac(1,1,1,1)=_E*_hardening_modulus/(_E+_hardening_modulus);
    }

    Stress.SetToZeros();
    Stress(1,1)=_TrialStress-_E*_DeltaGamma*Sign(_TrialStress);

    Mate.Rank2Materials("plastic_strain").SetToZeros();
    Mate.Rank2Materials("plastic_strain")(1,1)=MateOld.Rank2Materials("plastic_strain")(1,1)+_DeltaGamma*Sign(_TrialStress);

    Mate.ScalarMaterials("effective_plastic_strain")=MateOld.ScalarMaterials("effective_plastic_strain")+_DeltaGamma;


}

void Plastic1DMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) {

    //*********************************************************
    //*** get rid of unused warnings
    //*********************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}
 
    
    if(InputParams.size()<3){ 
        MessagePrinter::PrintErrorTxt("for 1d plastic material, three parameters are required, you need to give: E, yield stress, and hardening modulus");
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

