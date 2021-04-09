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
#pragma once

#include "MateSystem/MechanicsMaterialBase.h"

class LinearElasticMaterial: public MechanicsMaterialBase{
public:
    virtual void InitMaterialProperties(const double &t, const double &dt,const int &nDim,
                                        const Vector3d &gpCoord,const vector<double> &InputParams,
                                        const vector<double> &gpU,const vector<double> &gpUdot,
                                        const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUdot,
                                        Materials &Mate) override;

    virtual void ComputeMaterialProperties(const double &t, const double &dt,const int &nDim,
                                           const Vector3d &gpCoord,const vector<double> &InputParams,
                                           const vector<double> &gpU,const vector<double> &gpUOld,
                                           const vector<double> &gpUdot,const vector<double> &gpUdotOld,
                                           const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUOld,
                                           const vector<Vector3d> &gpGradUdot,const vector<Vector3d> &gpGradUdotOld,
                                           const Materials &MateOld, Materials &Mate) override;

private:
    virtual void ComputeStrain(const int &nDim,const vector<Vector3d> &GradDisp, RankTwoTensor &Strain) override;
    virtual void ComputeStressAndJacobian(const vector<double> &InputParams,const RankTwoTensor &Strain,RankTwoTensor &Stress,RankFourTensor &Jacobian) override;

private:
    RankTwoTensor _GradU,_Strain,_Stress,_I,_devStress;
    RankFourTensor _Jac;
};