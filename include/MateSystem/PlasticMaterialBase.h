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
//+++ Purpose: Define the base material for elasto-plastic materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include "MateSystem/BulkMaterialBase.h"

class PlasticMaterialBase: public BulkMaterialBase{
public:
    virtual void ComputeStrain(const int &nDim,const vector<Vector3d> &GradDisp,RankTwoTensor &Strain)=0;

    virtual double ComputeYieldFunction(const vector<double> &InputParams,const double &trial_strss,const double &effect_plastic_strain)=0;


};