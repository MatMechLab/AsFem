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
//+++ Purpose: implement the residual and jacobian for Miehe's
//+++          phase field fracture model
//+++          1) eta*dD/dt=2(1-D)H-(Gc/L)D+Gc*L*Lap(D) (see Eq. 47)
//+++          2) div(\Sigma)=0
//+++ Reference: A phase field model for rate-independent crack propagation:
//+++            Robust algorithmic implementation based on operator splits
//+++ DOI    : https://doi.org/10.1016/j.cma.2010.04.011
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ElmtSystem/BulkElmtBase.h"

/**
 * This class implement the calculation of Miehe's phase field fracture model
 */
class MieheFractureElmt:public BulkElmtBase{
public:
    /**
     * Function for different calculate action
     */
    virtual void ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
            const LocalElmtSolution &soln,const LocalShapeFun &shp,
            const Materials &Mate,const Materials &MateOld,
            ScalarMateType &gpProj,
            MatrixXd &localK,VectorXd &localR) override;

private:
    /**
     * This function calculate the residual of Miehe's phase field fracture model. <br>
     * \f$R_{u_{i}}^{I}=\sigma_{ij}N_{,j}^{I}\f$ <br>
     * \f$R_{d}=\eta\dot{d}N^{I}+2(d-1)\mathbb{H}N^{I}+\frac{G_{c}}{L}dN^{I}+G_{c}L\nabla d\nabla N^{I}\f$
     */
    virtual void ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) override;

    /**
     * This function calculate the jacobian matrix of Miehe's phase field fracture model. <br>
     * \f$K_{u_{i}u_{k}}^{IJ}=\frac{\partial R_{u_{i}}^{I}}{\partial u_{k}^{J}}=\mathbb{C}_{ijkl}N_{,j}^{I}N_{,l}^{J}\f$ <br>
     * \f$K_{dd}^{IJ}=\frac{\partial R_{d}^{I}}{\partial d^{J}}=\eta N^{J}N^{I}\mathrm{ctan[1]}+2N^{J}\mathcal{H}N^{I}+\frac{G_{c}}{L}N^{J}N^{I}+G_{c}L\nabla N^{J}\nabla N^{I}\f$
     */
    virtual void ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) override;

    /**
     * This function calculate the projected scalar variable for Miehe's phase field fracture model
     */
    virtual void ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) override;

};
