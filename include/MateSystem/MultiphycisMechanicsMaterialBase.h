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
//+++ Date   : 2021.11.14
//+++ Purpose: Define the base material for multiphysics coupled solid 
//+++          mechanics problem, all the coupled mechanical materials 
//+++          calculation should inherit this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include "MateSystem/BulkMaterialBase.h"

/**
 * The base class for solid mechanics materials, where the calculation of the stress \f$\mathbf{\sigma}\f$ and jacobian matrix \f$\mathbb{C}\f$ are required
 */
class MultiphysicsMechanicsMaterialBase: public BulkMaterialBase{
public:
    /**
     * Compute the deformation gradient tensor, it could be small strain \f$\mathbf{\varepsilon}\f$, finite strain deformation gradient tensor \f$\mathbf{F}=\mathbf{U}+\mathbf{I}\f$.
     * @param elmtinfo the current element information data
     * @param elmtsoln the curent element's solution, include the displacement and its gradient
     * @param F the rank-2 strain tensor
     */ 
    virtual void ComputeDeformationGradientTensor(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &F)=0;

    /**
     * Compute the stress \f$\mathbf{\sigma}\f$ and jacobian matrix \f$\mathbb{C}\f$ .
     * @param InputParams the input material parameters from your input file.
     * @param soln the local solution structure.
     * @param F the total deformation gradient tensor.
     * @param Strain the strain tensor.
     * @param Stress the calculated stress tensor \f$\mathbf{\sigma}\f$ for MechanicsElmt.
     * @param Jacobian the calculated 'elasticity' tensor \f$\mathbb{C}\f$, it is the rank-4 tensor for the general constitutive law.
     */
    virtual void ComputeConstitutiveLaws(const vector<double> &InputParams,const LocalElmtSolution &soln,const RankTwoTensor &F,RankTwoTensor &Strain,RankTwoTensor &Stress,RankFourTensor &Jacobian)=0;


};
