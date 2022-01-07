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
//+++ Purpose: calculate the stress and strain based on incremental
//+++          description for small strain case
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/IncrementSmallStrainMaterial.h"

void IncrementSmallStrainMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    // Here we do not consider any initial internal strains, stress
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}

    Mate.Rank2Materials("stress").SetToZeros();// we need the stress and strain from previous step
    Mate.Rank2Materials("strain").SetToZeros();
}

void IncrementSmallStrainMaterial::ComputeStrain(const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, RankTwoTensor &Strain){
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
//***********************************************************************************
void IncrementSmallStrainMaterial::ComputeStressAndJacobian(const vector<double> &InputParams,
                                                            const RankTwoTensor &Strain, RankTwoTensor &Stress,
                                                            RankFourTensor &Jacobian) {
    const double E=InputParams[0];
    const double nu=InputParams[1];

    Jacobian.SetToZeros();
    Jacobian.SetFromEandNu(E,nu);
    Stress=Jacobian.DoubleDot(Strain);
}
//*********************************************************************
void IncrementSmallStrainMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) {
    //***********************************************************************
    //*** get rid of unused warning
    //***********************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}

    _StrainOld=MateOld.Rank2Materials("strain");
    _StressOld=MateOld.Rank2Materials("stress");
    ComputeStrain(elmtinfo,elmtsoln,_Strain);
    _DeltaStrain=_Strain-_StrainOld;
    ComputeStressAndJacobian(InputParams,_DeltaStrain,_DeltaStress,_Jac);

    _Stress=_StressOld+_DeltaStress;

    _I.SetToZeros();
    _devStress=_Stress-_I*(_Stress.Trace()/3.0);

    Mate.ScalarMaterials("vonMises")=sqrt(1.5*_devStress.DoubleDot(_devStress));
    Mate.Rank2Materials("stress")=_Stress;
    Mate.Rank2Materials("strain")=_Strain;
    Mate.Rank4Materials("jacobian")=_Jac;

}
