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
//+++ Date   : 2021.04.09
//+++ Purpose: Calculate the material properties required by Miehe's
//+++          phase field fracture model
//+++           1) viscosity
//+++           2) Gc
//+++           3) L
//+++           4) H
//+++           5) dHdstrain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/PhaseFieldFractureMaterialBase.h"

class MieheFractureMaterial:public PhaseFieldFractureMaterialBase{
public:
    virtual void InitMaterialProperties(const int &nDim,const Vector3d &gpCoord,const vector<double> &InputParams,
                                        const vector<double> &gpU,const vector<double> &gpUdot,
                                        const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUdot,
                                        Materials &Mate) override;

    virtual void ComputeMaterialProperties(const double &t, const double &dt,const int &nDim, const Vector3d &gpCoord,
                                           const vector<double> &InputParams,
                                           const vector<double> &gpU,const vector<double> &gpUOld,
                                           const vector<double> &gpUdot,const vector<double> &gpUdotOld,
                                           const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUOld,
                                           const vector<Vector3d> &gpGradUdot,const vector<Vector3d> &gpGradUdotOld,
                                           const Materials &MateOld, Materials &Mate) override;

private:
    virtual void ComputeStrain(const int &nDim, const vector<Vector3d> &GradDisp,RankTwoTensor &strain) override;

    virtual void ComputeConstitutiveLaws(const vector<double> &InputParams,const RankTwoTensor &strain,const double &damage,
                                         const Materials &MateOld, Materials &Mate) override;

    virtual double DegradationFun(const double &x) override;
    virtual double DegradationFunDeriv(const double &x) override;

private:
    RankTwoTensor I,StressPos,StressNeg,DevStress,GradU;
    RankTwoTensor Strain,EpsPos,EpsNeg;
    RankFourTensor I4Sym,ProjPos,ProjNeg;

};