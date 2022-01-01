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
//+++ Date   : 2021.12.25
//+++ Purpose: implement the residual and jacobian for allen-cahh type
//+++          phase field fracture model
//+++          1) dd/dt=-M*(delta f/delta d)
//+++          2) div(\Sigma)=0
//+++ Reference: A continuum phase field model for fracture
//+++ DOI      : https://doi.org/10.1016/j.engfracmech.2010.08.009
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ElmtSystem/BulkElmtBase.h"

/**
 * This class implement the calculation of allen-cahn type phase field fracture model
 */
class AllenCahnFractureElmt:public BulkElmtBase{
public:
    /**
     * Function for different calculate action
     */
    virtual void ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
            const LocalElmtSolution &soln,const LocalShapeFun &shp,
            const Materials &Mate,const Materials &MateOld,
            ScalarMateType &gpProj,
            MatrixXd &localK,VectorXd &localR) override;

private:
    /**
     * This function calculate the residual of Miehe's phase field fracture model. <br>
     * \f$R_{u_{i}}^{I}=\sigma_{ij}N_{,j}^{I}\f$ <br>
     * \f$R_{d}=\eta\dot{d}N^{I}+M\frac{\delta\psi}{\delta d}N^{I}$
     */
    virtual void ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) override;

    /**
     * This function calculate the jacobian matrix of Miehe's phase field fracture model. <br>
     * \f$K_{u_{i}u_{k}}^{IJ}=\frac{\partial R_{u_{i}}^{I}}{\partial u_{k}^{J}}=\mathbb{C}_{ijkl}N_{,j}^{I}N_{,l}^{J}\f$ <br>
     * \f$K_{dd}^{IJ}=\frac{\partial R_{d}^{I}}{\partial d^{J}}=\eta N^{J}N^{I}\mathrm{ctan[1]}+M\frac{\delta^{2}\psi}{\delta d^{2}}N^{J}N^{I}\f$
     */
    virtual void ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) override;

    /**
     * This function calculate the projected scalar variable for allen-cahn type phase field fracture model
     */
    virtual void ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) override;

};
