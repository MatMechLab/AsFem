//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
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

/**
 * This class defines the base class for phase-field fracture model
 */
class PhaseFieldFractureMaterialBase: public BulkMaterialBase{
protected:
    /**
     * Compute the strain, it could be small strain \f$\mathbf{\varepsilon}\f$, Green-Lagrange tensor \f$\mathbf{E}=\frac{1}{2}(\mathbf{F}^{T}\mathbf{F}-\mathbf{I})\f$.
     * @param elmtinfo the current element information data
     * @param elmtsoln the curent element's solution, include the displacement and its gradient
     * @param Strain the rank-2 strain tensor
     */ 
    virtual void ComputeStrain(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &Strain)=0;

    /**
     * Compute the stress \f$\mathbf{\sigma}\f$ and jacobian matrix \f$\mathbb{C}\f$ .
     * @param InputParams the input material parameters from your input file.
     * @param damage the order parameter for damage field.
     * @param Strain the input strain tensor, it is calculated from ComputeStrain function.
     * @param Stress the calculated stress tensor \f$\mathbf{\sigma}\f$ for MechanicsElmt.
     * @param Jacobian the calculated 'elasticity' tensor \f$\mathbb{C}\f$, it is the rank-4 tensor for the general constitutive law.
     */
    virtual void ComputeConstitutiveLaws(const vector<double> &InputParams,const double &damage,const RankTwoTensor &strain,RankTwoTensor &Stress,RankFourTensor &Jacobian,const Materials &MateOld, Materials &Mate)=0;
    
    /**
     * The degradation function
     * @param x the damage variable value
     */ 
    virtual double DegradationFun(const double &x)=0;
    
    /**
     * The degradation function
     * @param x the damage variable value
     */
    virtual double DegradationFunDeriv(const double &x)=0;

};
