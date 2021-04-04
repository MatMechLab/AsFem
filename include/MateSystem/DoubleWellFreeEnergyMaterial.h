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
//+++ Purpose: Calculate the free energy, chemical potential and its
//+++          derivatives of double well free energy material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/FreeEnergyMaterialBase.h"

class DoubleWellFreeEnergyMaterial: public FreeEnergyMaterialBase{
public:
    DoubleWellFreeEnergyMaterial();
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

protected:
    virtual void ComputeF(const vector<double> &InputParams,const vector<double> &U,const vector<double> &dUdt,vector<double> &F) override;
    virtual void ComputedFdU(const vector<double> &InputParams,const vector<double> &U,const vector<double> &dUdt, vector<double> &dF) override;
    virtual void Computed2FdU2(const vector<double> &InputParams,const vector<double> &U,const vector<double> &dUdt, vector<double> &d2F) override;

private:
    double c;
    vector<double> _F,_dFdc,_d2Fdc2;
};

