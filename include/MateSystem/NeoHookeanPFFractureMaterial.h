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
//+++ Date   : 2022.10.28
//+++ Purpose: Calculate the material properties defined in neohookean
//+++          type finite strain phase-field fracture model
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/BulkMaterialBase.h"
#include "MateSystem/FreeEnergyMaterialBase.h"
#include "MateSystem/ElasticMaterialBase.h"

/**
 * This class calculate the material properties for the allen-cahn fracture element
 */
class NeoHookeanPFFractureMaterial:public BulkMaterialBase,
                                   public FreeEnergyMaterialBase,
                                   public ElasticMaterialBase{
public:
    NeoHookeanPFFractureMaterial();
    ~NeoHookeanPFFractureMaterial();
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
    /**
     * Calculate the free energy and its first/second order (partial) derivatives.
     * @param parameters the json parameters read from input file
     * @param args the variables for the free energy expression, it could be concentration or order parameters
     * @param F the system free energy
     * @param dFdargs the first order derivatives of F with respect to its own args
     * @param d2Fdargs2 the second order (partial) derivatives of F with respect to different args, off-diagnoal part for partial derivatives
    */
    virtual void computeFreeEnergyAndDerivatives(const nlohmann::json &parameters,
                                                 const VectorXd &args,
                                                 VectorXd       &F,
                                                 VectorXd       &dFdargs,
                                                 MatrixXd       &d2Fdargs2) override;

private:
    /**
     * The degradation function
     * @param x the damage variable, 0-> undamage,1->full damage
    */
    double g(const double &x){
        return (1.0-x)*(1.0-x);
    }

    /**
     * The 1st order derivative of the degradation function
     * @param x the damage variable, 0-> undamage,1->full damage
    */
    double dg(const double &x){
        return 2.0*(x-1.0);
    }

    /**
     * The 2nd order derivative of the degradation function
     * @param x the damage variable, 0-> undamage,1->full damage
    */
    double d2g(const double &x){
        if(x){}
        return 2.0;
    }

private:
    double m_d;/**< the damage variable */
    double m_psi,m_psipos,m_psineg;/**< the different free energies */
    double m_stabilizer=1.0e-5;/**< stabilizer to get rid of rigid body motion in fully damaged region */

    VectorXd m_args;/**< the variables */
    VectorXd m_F;/**< the system free energy*/
    VectorXd m_dFdargs;/**< the first order derivatives of F */
    MatrixXd m_d2Fdargs2;/**< the second order derivatives of F */

    double m_Je,m_Je23;/**< determinte of deformation tensor */
    double m_I1,m_I1bar;/**< first invariant*/
    Rank2Tensor m_GradU;/** the displacement's gradient */
    Rank2Tensor m_Estrain;/**< the Lagrange-Green strain tensor */
    Rank2Tensor m_Fe;/**< the mechanical deformation tensor */
    Rank2Tensor m_Ce,m_CeInv;/**< the right Cauchy-Green tensor and its inverse */
    Rank2Tensor m_PK2stress,m_PK2stress_pos,m_PK2stress_neg;/**< the pk2 stress and its positive and negative part */
    Rank2Tensor m_dPK2stress_dD;/**< the derivative of pk2 stress w.r.t. damage */

    Rank2Tensor m_PK1stress,m_dPK1stress_dD;/**< the PK1 stress and its derivatives */
    Rank2Tensor m_devstress;/**< the deviatoric part of the pk1 stress */

    Rank2Tensor m_I;/**< the identity tensor */

    Rank4Tensor m_jacobian,m_jacobian_pos,m_jacobian_neg;/**< the jacobian tensor and its positive and negative part */

};