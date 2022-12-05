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
//+++ Date   : 2022.08.22
//+++ Purpose: Implement the calculation of Saint Venant
//+++          hyperelastic material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/BulkMaterialBase.h"
#include "MateSystem/ElasticMaterialBase.h"


/**
 * This class implement the constitutive law for linear elastic materials
 */
class SaintVenantMaterial:public BulkMaterialBase,
                          public ElasticMaterialBase{
protected:
    /**
     * Initial the preset material properties, if you don't need the history information of some materials, then you can avoid calling this function
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
     * compute the elastic strain tensor
     * @param dim the dimension of current analysis
     * @param gradU the gradient of displacement
     * @param strain the output elastic strain tensor
     */
    virtual void computeStrain(const int &dim,const Rank2Tensor &gradU,Rank2Tensor &strain) override;
    /**
     * compute the stress and jacobian tensor
     * @param t_params the json parameters defined in the input file
     * @param dim the dimension of current analysis
     * @param strain the strain tensor
     * @param stress the output stress tensor
     * @param jacobian the output jacobian tensor
     */
    virtual void computeStressAndJacobian(const nlohmann::json &t_params,
                                          const int &dim,
                                          const Rank2Tensor &strain,
                                          Rank2Tensor &stress,
                                          Rank4Tensor &jacobian) override;

private:
    double m_Psi;/**< the elastic strain energy density */
    Rank2Tensor m_devStress,m_devStrain;/**< for the deviatoric part of stress and strain */
    Rank2Tensor m_gradU;/**< the displacement gradient */
    Rank2Tensor m_F,m_I;/**< for the deformation tensor and identity tensor */
    Rank2Tensor m_C;/**< for the right Cauchy-Green tensor */
    Rank2Tensor m_strain,m_pk2_stress,m_stress;/**< local stress and strain tensor */
    Rank4Tensor m_jacobian,m_I4Sym;/**< local jacobian tensor */

};