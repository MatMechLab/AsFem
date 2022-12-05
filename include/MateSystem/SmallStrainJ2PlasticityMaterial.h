//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.11.13
//+++ Purpose: Implement the J2 plasticity material in small strain case
//+++          the algorithm is taken from Prof. Simo's book
//+++          details can be found in BOX 3.1(P122), BOX 3.2(P124)
//+++ Ref    : Computational Inelasticity, Book by J. C. Simo and Thomas J.R. Hughes
//+++ DOI    : https://link.springer.com/book/10.1007/b98904
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/BulkMaterialBase.h"
#include "MateSystem/PlasticMaterialBase.h"


/**
 * This class calculate the constitituve laws for J2 plasticity material in small strain case.
*/
class SmallStrainJ2PlasticityMaterial:public BulkMaterialBase,
                                      public PlasticMaterialBase{
public:
    SmallStrainJ2PlasticityMaterial();
    ~SmallStrainJ2PlasticityMaterial();
protected:
    /**
     * Initial the preset material properties, if you don't need the history information of some materials, 
     * then you can avoid calling this function
     * @param t_inputparams the input material parameters read from the json file
     * @param t_elmtinfo the data structure for the local element information
     * @param t_elmtsoln the solutions, i.e., U and V of the local element
     * @param Mate the materials (container) to be initialized
     */
    virtual void initMaterialProperties(const nlohmann::json &t_inputparams,
                                        const LocalElmtInfo &t_elmtinfo,
                                        const LocalElmtSolution &t_elmtsoln,
                                        MaterialsContainer &t_mate) override;

    /**
     * Compute the material property accroding to your model
     * @param t_inputparams the input material parameters read from the input file
     * @param t_elmtinfo the data structure for the local element information
     * @param t_elmtsoln the solutions, i.e., U and V of the local element
     * @param t_mateold the materials from previous step
     * @param t_mate the materials to be calculated
     */
    virtual void computeMaterialProperties(const nlohmann::json &t_inputparams,
                                           const LocalElmtInfo &t_elmtinfo,
                                           const LocalElmtSolution &t_elmtsoln,
                                           const MaterialsContainer &t_mateold,
                                           MaterialsContainer &t_mate) override;

private:
    /**
     * Evaluate the value of the yield function
     * @param t_parameters the input json parameters read from input file
     * @param t_args the input variables, i.e., effective plastic strain, etc.
     * @param t_mate the input material, i.e. cauchy stress, Mandel stress, plastic strain, etc.
    */
    virtual double computeYieldFunction(const nlohmann::json &t_parameters,
                                        const VectorXd &t_args,
                                        const MaterialsContainer &t_mate) override;
    /**
     * Evaluate the 1st order derivative value of the yield function
     * @param t_parameters the input json parameters read from input file
     * @param t_args the input variables, i.e., effective plastic strain, etc.
     * @param t_mate the input material, i.e. cauchy stress, Mandel stress, plastic strain, etc.
    */
    virtual double computeYieldFunctionDeriv(const nlohmann::json &t_parameters,
                                             const VectorXd &t_args,
                                             const MaterialsContainer &t_mate) override;
    
    /**
     * Evalute the admissible stresses status during the elasto-plastic deformation, 
     * where one can use the radial return mapping for the calculation.
     * @param t_parameters the input json parameters read from input file
     * @param t_elmtinfo the local element info structure
     * @param t_elmtsoln the local element solution 
     * @param t_mateold the old material constain from previous step
     * @param t_total_strain the total strain or total deformation gradient tensor
     * @param t_mate the current material container
    */
    virtual void computeAdmissibleStressState(const nlohmann::json &t_parameters,
                                              const LocalElmtInfo &t_elmtinfo,
                                              const LocalElmtSolution &t_elmtsoln,
                                              const MaterialsContainer &t_mateold,
                                              const Rank2Tensor &t_total_strain,
                                              MaterialsContainer &t_mate) override;


private:
    double m_eff_plastic_strain;/**< the effective plastic strain */
    double m_F;/**< Yield function */
    double m_gamma,m_theta,m_theta_bar;
    Rank2Tensor m_total_strain,m_GradU;/** total strain tensor and disp gradient tensor */
    Rank2Tensor m_plastic_strain_old;/**< the plastic strain of previous step */
    Rank2Tensor m_stress_trial,m_N;/**< the trial stress and plastic strain tensor */
    Rank2Tensor m_devStress,m_stress,m_dev_strain,m_I;
    Rank4Tensor m_I4Sym;/**< symmetric identity-4 tensor */

private:
    VectorXd m_args;/**< the arguments for different dofs */

};