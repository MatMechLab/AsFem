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
//+++ Purpose: implement the residual and jacobian for Poisson equation
//+++          Sigma*div(grad(phi))=F
//+++          Both Sigma and F could be nonlinear function !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ElmtSystem/BulkElmtBase.h"


/**
 * This class implement the calculation of poisson equation
 */
class PoissonElmt:public BulkElmtBase{
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
     * This function calculate the residual of the poisson equation. <br>
     * \f$R_{\phi}^{I}=\nabla\phi\nabla N^{I}+FN^{I}\f$
     */
    virtual void ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) override;

    /**
     * This function calculate the jacobian matrix of the poisson equation. <br>
     * \f$K_{\phi\phi}^{IJ}=\sigma\nabla N^{J}\nabla N^{I}\f$
     */
    virtual void ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) override;
    
    /**
     * This function calculate the projected scalar variable for the poisson equation
     */
    virtual void ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) override;

};
