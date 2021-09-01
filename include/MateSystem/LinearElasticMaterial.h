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

/**
 * This class implement the isotropic linear elastic material for the small strain case.
 */
class LinearElasticMaterial: public MechanicsMaterialBase{
public:
    /**
     * Initialize material properties in LinearElasticMaterial
     */
    virtual void InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate) override;
   
    /**
     * Compute the stress and jacobian in LinearElasticMaterial
     */ 
    virtual void ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) override;

private:
    /**
     * compute small strain
     */
    virtual void ComputeStrain(const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, RankTwoTensor &Strain) override;

    /**
     * compute the stress and jacobian for small deformation case.
     */
    virtual void ComputeStressAndJacobian(const vector<double> &InputParams, const RankTwoTensor &Strain, RankTwoTensor &Stress, RankFourTensor &Jacobian) override;


private:
    RankTwoTensor _GradU,_Strain,_Stress,_I,_devStress;
    RankFourTensor _Jac;
};
