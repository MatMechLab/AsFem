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
//+++ Purpose: define the base class for phase field fracture model
//+++          related material properties calculation
//+++          All the phase field fracture material should inherit
//+++          this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/BulkMaterialBase.h"
#include "Utils/MathFuns.h"

class PhaseFieldFractureMaterialBase: public BulkMaterialBase{
protected:
    virtual void ComputeStrain(const int &nDim,const vector<Vector3d> &GradDisp,RankTwoTensor &Strain)=0;
    virtual void ComputeConstitutiveLaws(const vector<double> &InputParams,const RankTwoTensor &strain,const double &damage,
                                         const Materials &MateOld, Materials &Mate)=0;
    virtual double DegradationFun(const double &x)=0;
    virtual double DegradationFunDeriv(const double &x)=0;

};