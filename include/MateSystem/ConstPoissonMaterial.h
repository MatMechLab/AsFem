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
//+++ Date   : 2021.04.10
//+++ Purpose: Calculate the material properties required by Poisson
//+++          element. In this code, we can define:
//+++           1) Sigma
//+++           2) dSigma/du(=0)
//+++           3) F
//+++           4) dF/du(=0)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/BulkMaterialBase.h"

/**
 * This class calculate the constant sigma and F for the poisson equation
 * \f$\sigma\nabla^{2}\phi=F\f$
 */
class ConstPoissonMaterial:public BulkMaterialBase{
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
