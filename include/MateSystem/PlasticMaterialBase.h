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

/**
 * This class define the general plasticity material calculation
 */
class PlasticMaterialBase: public BulkMaterialBase{
public:
    /**
     * Compute the strain, it could be small strain \f$\mathbf{\varepsilon}\f$, Green-Lagrange tensor \f$\mathbf{E}=\frac{1}{2}(\mathbf{F}^{T}\mathbf{F}-\mathbf{I})\f$.
     * @param elmtinfo the current element information data
     * @param elmtsoln the curent element's solution, include the displacement and its gradient
     * @param Strain the rank-2 strain tensor
     */ 
    virtual void ComputeStrain(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &Strain)=0;

    /**
     * Compute the yield function for elasto-plastic check
     * @param InputParams the material parameters read from the input file
     * @trial_strss the trial stress for plastic check(the yield function)
     * @effect_plastic_strain_old the previous effective plastic strain
     */
    virtual double ComputeYieldFunction(const vector<double> &InputParams,const double &trial_strss,const double &effect_plastic_strain_old)=0;

    /**
     * Compute the admissible stress state
     * @param InputParams the material parameters read from input file 
     * @param elmtinfo the current element's local information
     * @param elmtsoln the solution of current element, i.e., disp and its gradients
     * @param Strain the strain tensor
     * @param MateOld the materials from previous step
     * @param Mate the materials of current step(to be calculated)
     * @param Stress computed stress
     * @param Jac the rank-4 tensor for the constitutive law
     */
    virtual void ComputeAdmissibleStressState(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,const RankTwoTensor &Strain,const Materials &MateOld,Materials &Mate,RankTwoTensor &Stress,RankFourTensor &Jac)=0;


};
