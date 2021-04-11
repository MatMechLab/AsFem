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
//+++ Purpose: calculate the stress and strain based on incremental
//+++          description for small strain case
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/IncrementSmallStrainMaterial.h"

void IncrementSmallStrainMaterial::InitMaterialProperties(const int &nDim, const Vector3d &gpCoord,
                                                          const vector<double> &InputParams, const vector<double> &gpU,
                                                          const vector<double> &gpUdot, const vector<Vector3d> &gpGradU,
                                                          const vector<Vector3d> &gpGradUdot, Materials &Mate) {
    //*************************************************************************
    //*** get rid of unused warning
    //*************************************************************************
    if(nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUdot[0]||
       gpGradU[0](1)||gpGradUdot[0](1)||Mate.ScalarMaterials.size()){}

    Mate.Rank2Materials["stress"].SetToZeros();// we need the stress and strain from previous step
    Mate.Rank2Materials["strain"].SetToZeros();
}
//**************************************************************************************
void IncrementSmallStrainMaterial::ComputeStrain(const int &nDim, const vector<Vector3d> &GradDisp,
                                                 RankTwoTensor &Strain) {
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
void IncrementSmallStrainMaterial::ComputeMaterialProperties(const double &t, const double &dt, const int &nDim,
                                                             const Vector3d &gpCoord, const vector<double> &InputParams,
                                                             const vector<double> &gpU, const vector<double> &gpUOld,
                                                             const vector<double> &gpUdot,
                                                             const vector<double> &gpUdotOld,
                                                             const vector<Vector3d> &gpGradU,
                                                             const vector<Vector3d> &gpGradUOld,
                                                             const vector<Vector3d> &gpGradUdot,
                                                             const vector<Vector3d> &gpGradUdotOld,
                                                             const Materials &MateOld, Materials &Mate) {
    //***********************************************************************
    //*** get rid of unused warning
    //***********************************************************************
    if(t||dt||gpCoord(1)||gpU[0]||gpUOld[0]||
       gpUdot[0]||gpUdotOld[0]||gpGradU[0](1)||gpGradUOld[0](1)||
       gpGradUdot[0](1)||gpGradUdotOld[0](1)||MateOld.ScalarMaterials.size()){}

    _StrainOld=MateOld.Rank2Materials.at("strain");
    _StressOld=MateOld.Rank2Materials.at("stress");
    ComputeStrain(nDim,gpGradU,_Strain);
    _DeltaStrain=_Strain-_StrainOld;
    ComputeStressAndJacobian(InputParams,_DeltaStrain,_DeltaStress,_Jac);

    _Stress=_StressOld+_DeltaStress;

    _I.SetToZeros();
    _devStress=_Stress-_I*(_Stress.Trace()/3.0);

    Mate.ScalarMaterials["vonMises"]=sqrt(1.5*_devStress.DoubleDot(_devStress));
    Mate.Rank2Materials["stress"]=_Stress;
    Mate.Rank2Materials["strain"]=_Strain;
    Mate.Rank4Materials["jacobian"]=_Jac;
}