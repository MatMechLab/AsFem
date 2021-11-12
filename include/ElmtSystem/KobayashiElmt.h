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
//+++ Date   : 2021.11.12
//+++ Purpose: implement the residual and jacobian for Kobayashi's
//+++          dendrite model
//+++          the governing equations are:
//+++          1) deta/dt=L*div(k*dk*v)+L*div(k*k*\nabla\eta)-L*df/deta
//+++          2) dT/dt=Lap(T)+K*deta/dt             
//+++ Dofs   : 1. eta, 2, T
//+++ Ref    : https://doi.org/10.1016/0167-2789(93)90120-P
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ElmtSystem/BulkElmtBase.h"

/**
 * This class implement the calculation of Kobayashi's dendrite equations.
 * The implementation requires two Dofs, one for the order parameter, another one is the temperature.
 */
class KobayashiElmt:public BulkElmtBase{
    
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
     * This function is responsible for the residual calculation of Kobayashi equation.
     * The explanation of each arguments can be found in ElmtBase class.<br>
     */
    virtual void ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) override;

    /**
     * This function calculate the jacobian matrix for Kobayashi equation.<br>
     */
    virtual void ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) override;

    /**
     * For the projected scalar variables in Kobayashi element 
     */
    virtual void ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) override;

private:
    Vector3d V,dV;/**< used in the calculaton of residual and jacobian >*/
    double L,K,dK,Latent;
    double dFdeta,d2Fdeta2,d2FdetadT;
    Vector3d dKdGradEta,ddKdGradEta;


};
