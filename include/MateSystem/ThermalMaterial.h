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
//+++ Date   : 2022.01.29
//+++ Purpose: Calculate the material properties required by thermal
//+++          conduct element. In this code, we can define:
//+++           1) rho--->density
//+++           2) cp --->capacity
//+++           3) qv --->body heat source
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/BulkMaterialBase.h"

/**
 * This class calculate the thermal material properties for heat equation
 */
class ThermalMaterial:public BulkMaterialBase{
public:
    /**
     * Initialze the material properties
     */
    virtual void InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate) override;

    /**
     * Calculate the \f$\sigma\f$, \f$F\f$ and their derivative(=0 for constant case)
     */
    virtual void ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) override;

};
