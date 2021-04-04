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
//+++ Purpose: Define the base material for free energy materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include "MateSystem/BulkMaterialBase.h"

class FreeEnergyMaterialBase: public BulkMaterialBase{
public:
    virtual void ComputeF(const vector<double> &InputParams,const vector<double> &U,const vector<double> &dUdt,vector<double> &F)=0;
    virtual void ComputedFdU(const vector<double> &InputParams,const vector<double> &U,const vector<double> &dUdt,vector<double> &dF)=0;
    virtual void Computed2FdU2(const vector<double> &InputParams,const vector<double> &U,const vector<double> &dUdt,vector<double> &d2F)=0;

};