//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.04.10
//+++ Purpose: implement the residual and jacobian for cahn-hilliard 
//+++          equation
//+++          dc/dt=div(M*grad(u))
//+++             mu=deltaF/deltaC
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ElmtSystem/BulkElmtBase.h"


/**
 * This class implement the calculation of cahn-hilliard equation
 */
class CahnHilliardElement:public BulkElmtBase{
public:
    /**
     * This function responsible for calculation of the residual and jacobian,
     * @param CalcType the calculation type of FEM analysis, i.e., residual-calc, jacobian-calc, projection-calc
     * @param ElmtInfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param Ctan 1x3 vector, where ctan[0] is responsible for the non-time-derivative part in the K matrix, while ctan[1] and ctan[2] represent the coeffecient for the 1st and 2nd order time derivatives in the K matrix
     * @param ElmtSoln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param Shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param MateOld the materials of previous step
     * @param Mate the materials of current step
     * @param LocalK the K matrix of local element, it could be 1x1 or 5x5 depends on your local nDOfs
     * @param LocalR the Residual vector of local element,it could 1 or 5 depends on your local nDofs
     */
    virtual void computeAll(const FECalcType &CalcType,
                            const LocalElmtInfo &ElmtInfo,
                            const double (&Ctan)[3],
                            const LocalElmtSolution &ElmtSoln,
                            const LocalShapeFun &Shp,
                            const MaterialsContainer &MateOld,
                            const MaterialsContainer &Mate,
                            MatrixXd &LocalK,VectorXd &LocalR) override;

protected:
    /**
     * This function is responsible for the local residual vector calculation
     * @param ElmtInfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param ElmtSoln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param Shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param MateOld the materials of previous step
     * @param Mate the materials of current step
     * @param LocalR the Residual vector of local element,it could 1 or 5 depends on your local nDofs 
     */
    virtual void computeResidual(const LocalElmtInfo &ElmtInfo,
                                 const LocalElmtSolution &ElmtSoln,
                                 const LocalShapeFun &Shp,
                                 const MaterialsContainer &MateOld,
                                 const MaterialsContainer &Mate,
                                 VectorXd &LocalR) override;

    /**
     * This function responsible for the  jacobian matrix calculation
     * @param ElmtInfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param Ctan 1x3 vector, where ctan[0] is responsible for the non-time-derivative part in the K matrix, while ctan[1] represents the coeffecient for the time derivatives in the K matrix
     * @param ElmtSoln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param Shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param MateOld the materials of previous step
     * @param Mate the materials of current step
     * @param LocalK the K matrix of local element, it could be 1x1 or 5x5 depends on your local nDOfs
     */
    virtual void computeJacobian(const LocalElmtInfo &ElmtInfo,
                                 const double (&Ctan)[3],
                                 const LocalElmtSolution &ElmtSoln,
                                 const LocalShapeFun &Shp,
                                 const MaterialsContainer &MateOld,
                                 const MaterialsContainer &Mate,
                                 MatrixXd &LocalK) override;
};