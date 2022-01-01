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
//+++ Purpose: implement the residual and jacobian for general
//+++          Cahn-Hilliard equation
//+++          dc/dt=div(M*grad(mu))
//+++             mu=df/dc-Kappa*Delta(c)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ElmtSystem/BulkElmtBase.h"

/**
 * This class implement the calculation of CahnHilliard equation.
 * The implementation requires two Dofs, one for the concentration, another one is the chemical potential, the mixed form is used to solve this 4th order PDE.
 */
class CahnHilliardElmt:public BulkElmtBase{
    
public:
    /**
     * The function for different calc action
     */
    virtual void ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
            const LocalElmtSolution &soln,const LocalShapeFun &shp,
            const Materials &Mate,const Materials &MateOld,
            ScalarMateType &gpProj,
            MatrixXd &localK,VectorXd &localR) override;

private:

    /**
     * This function is responsible for the residual calculation of CahnHilliard equation.
     * The explanation of each arguments can be found in ElmtBase class.<br>
     * \f$R_{c}^{I}=\dot{c}N^{I}+M\nabla c\nabla\mu\f$ <br>
     * \f$R_{\mu}^{I}=\mu N^{I}-\frac{\partial\psi}{\partial c}N^{I}-\kappa\Delta c\f$
     */
    virtual void ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) override;

    /**
     * This function calculate the jacobian matrix for CahnHilliard equation.<br>
     * \f$K_{cc}^{IJ}=\frac{\partial R_{c}^{I}}{\partial c^{J}}\f$ <br>
     * \f$K_{c\mu}^{IJ}=\frac{\partial R_{c}^{I}}{\partial\mu^{J}}\f$ <br>
     * \f$K_{\mu c}^{IJ}=\frac{\partial R_{\mu}^{I}}{\partial c^{J}}\f$ <br>
     * \f$K_{\mu\mu}^{IJ}=\frac{\partial R_{\mu}^{I}}{\partial\mu^{J}}\f$ <br>
     */
    virtual void ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) override;

    /**
     * For the projected scalar variables in CahnHilliard element 
     */
    virtual void ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) override;

};
