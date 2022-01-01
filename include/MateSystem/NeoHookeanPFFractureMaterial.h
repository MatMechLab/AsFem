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
//+++ Date   : 2021.11.07
//+++ Purpose: Calculate the material properties required by Miehe's
//+++          phase field fracture model by using neohookean material
//+++           1) viscosity
//+++           2) Gc
//+++           3) L
//+++           4) H
//+++           5) dHdstrain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/PhaseFieldFractureMaterialBase.h"

/**
 * This class implement the calculation for the constitutive laws based Miehe's phase field fracture model, where the constitutive laws are calculated for neohookean material.
 */
class NeoHookeanPFFractureMaterial:public PhaseFieldFractureMaterialBase{
public:
    /**
     * Initialize material properties in Miehe's phase field fracture model
     */
    virtual void InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate) override;
   
    /**
     * Compute the stress and jacobian in Miehe's phase field fracture model
     */ 
    virtual void ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) override;

private:

    /**
     * Compute the strain, it could be small strain \f$\mathbf{\varepsilon}\f$, Green-Lagrange tensor \f$\mathbf{E}=\frac{1}{2}(\mathbf{F}^{T}\mathbf{F}-\mathbf{I})\f$.
     * @param elmtinfo the current element information data
     * @param elmtsoln the curent element's solution, include the displacement and its gradient
     * @param Strain the rank-2 strain tensor
     */ 
    virtual void ComputeStrain(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &Strain) override;

    /**
     * Compute the stress \f$\mathbf{\sigma}\f$ and jacobian matrix \f$\mathbb{C}\f$ .
     * @param InputParams the input material parameters from your input file.
     * @param damage the order parameter for damage field.
     * @param Strain the input strain tensor, it is calculated from ComputeStrain function.
     * @param Stress the calculated stress tensor \f$\mathbf{\sigma}\f$ for MechanicsElmt.
     * @param Jacobian the calculated 'elasticity' tensor \f$\mathbb{C}\f$, it is the rank-4 tensor for the general constitutive law.
     */
    virtual void ComputeConstitutiveLaws(const vector<double> &InputParams,const double &damage,const RankTwoTensor &strain,RankTwoTensor &Stress,RankFourTensor &Jacobian,const Materials &MateOld, Materials &Mate) override;
    
    /**
     * The degradation function
     * @param x the damage variable value
     */ 
    virtual double DegradationFun(const double &x) override;
    
    /**
     * The degradation function
     * @param x the damage variable value
     */
    virtual double DegradationFunDeriv(const double &x) override;

private:
    RankTwoTensor _I,_Stress,_PK2,_CauchyStress,_DevStress,_GradU;
    RankTwoTensor _Strain,_F,_C,_Cinv,_E;
    RankFourTensor _Jacobian;
    RankTwoTensor _StressPos,_StressNeg;
    RankFourTensor _JacPos,_JacNeg,I4;
    double _J,_J23,_I1,_I1bar;

};
