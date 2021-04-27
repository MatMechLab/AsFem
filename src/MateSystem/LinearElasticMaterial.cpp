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
//+++ Date   : 2021.04.04
//+++ Purpose: Implement the calculation for linear elastic material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/LinearElasticMaterial.h"

void LinearElasticMaterial::InitMaterialProperties(const int &nDim,const Vector3d &gpCoord,const vector<double> &InputParams,
                                                   const vector<double> &gpU,const vector<double> &gpUdot,
                                                   const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUdot,
                                                   Materials &Mate) {
    // Here we do not consider any initial internal strains, stress
    if(nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUdot[0]||
       gpGradU[0](1)||gpGradUdot[0](1)||Mate.ScalarMaterials.size()){

    }
}


void LinearElasticMaterial::ComputeStrain(const int &nDim,const vector<Vector3d> &GradDisp, RankTwoTensor &Strain) {
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

void LinearElasticMaterial::ComputeStressAndJacobian(const vector<double> &InputParams,const RankTwoTensor &Strain,
                                                     RankTwoTensor &Stress,RankFourTensor &Jacobian) {

    const double E=InputParams[0];
    const double nu=InputParams[1];

    Jacobian.SetToZeros();
    Jacobian.SetFromEandNu(E,nu);
    Stress=Jacobian.DoubleDot(Strain);
}

void LinearElasticMaterial::ComputeMaterialProperties(const double &t, const double &dt,const int &nDim,
                                                      const Vector3d &gpCoord,const vector<double> &InputParams,
                                                      const vector<double> &gpU,const vector<double> &gpUOld,
                                                      const vector<double> &gpUdot,const vector<double> &gpUdotOld,
                                                      const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUOld,
                                                      const vector<Vector3d> &gpGradUdot,const vector<Vector3d> &gpGradUdotOld,
                                                      const Materials &MateOld, Materials &Mate) {

    if(t||dt||gpCoord(1)||gpU[0]||gpUOld[0]||
    gpUdot[0]||gpUdotOld[0]||gpGradU[0](1)||gpGradUOld[0](1)||
    gpGradUdot[0](1)||gpGradUdotOld[0](1)||MateOld.ScalarMaterials.size()){}// get rid of unused warning

    if(InputParams.size()<2){
        MessagePrinter::PrintErrorTxt("for linear elastic material, two parameters are required, you need to give: E and nu");
        MessagePrinter::AsFem_Exit();
    }

    ComputeStrain(nDim,gpGradU,_Strain);
    ComputeStressAndJacobian(InputParams,_Strain,_Stress,_Jac);

    _I.SetToIdentity();
    _devStress=_Stress-_I*(_Stress.Trace()/3.0);
    Mate.ScalarMaterials["vonMises"]=sqrt(1.5*_devStress.DoubleDot(_devStress));
    Mate.Rank2Materials["strain"]=_Strain;
    Mate.Rank2Materials["stress"]=_Stress;
    Mate.Rank4Materials["jacobian"]=_Jac;

}
