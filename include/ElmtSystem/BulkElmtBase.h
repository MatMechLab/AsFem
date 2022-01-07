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
//+++ Date   : 2021.04.09
//+++ Purpose: here we define the base class for the FEM calculation
//+++          i.e., the residual and jacobian calculation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>

#include "Utils/Vector3d.h"
#include "Utils/MathFuns.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"

#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

#include "FESystem/FECalcType.h"
#include "MateSystem/Materials.h"
#include "ElmtSystem/LocalElmtData.h"

#include "Utils/MessagePrinter.h"

/**
 * This abstract class defines the basic calculation of the local elmenent.
 * It should be noted that, this calculation is focusing on single 
 * quadrature point.<br>
 * For example, the residual looks like:\f$\int_{\Omega}\nabla c\nabla N^{I}dV\f$, the bulk 
 * element only calculate the \f$\nabla c\nabla N^{I}\f$ part.
 */
class BulkElmtBase{
public:
    /**
     * This function responsible for the residual, jacobian, and projection variable calculation, 
     * all the child class,i.e., mechanics and cahnhilliard, should override this function explicitly!!!
     * @param calctype the calculation type of FEM analysis, i.e., residual-calc, jacobian-calc, projection-calc
     * @param elmtinfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param ctan 1x2 vector, where ctan[0] is responsible for the non-time-derivative part in the K matrix, while ctan[1] represents the coeffecient for the time derivatives in the K matrix
     * @param soln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param Mate the materials of current step
     * @param MateOld the materials of previous step
     * @param gpProj the projected scalar variable in current element
     * @param localK the K matrix of local element, it could be 1x1 or 5x5 depends on your local nDOfs
     * @param localR the Residual vector of local element,it could 1 or 5 depends on your local nDofs
     */
    virtual void ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
            const LocalElmtSolution &soln,const LocalShapeFun &shp,
            const Materials &Mate,const Materials &MateOld,
            ScalarMateType &gpProj,
            MatrixXd &localK,VectorXd &localR)=0;

protected:
    /**
     * This function is responsible for the local residual vector calculation
     * @param elmtinfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param soln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param Mate the materials of current step
     * @param MateOld the materials of previous step
     * @param localR the Residual vector of local element,it could 1 or 5 depends on your local nDofs 
     */
    virtual void ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR)=0;

    /**
     * This function responsible for the  jacobian matrix calculation
     * @param elmtinfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param ctan 1x2 vector, where ctan[0] is responsible for the non-time-derivative part in the K matrix, while ctan[1] represents the coeffecient for the time derivatives in the K matrix
     * @param soln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param Mate the materials of current step
     * @param MateOld the materials of previous step
     * @param gpProj the projected scalar variable in current element
     * @param localK the K matrix of local element, it could be 1x1 or 5x5 depends on your local nDOfs
     */
    virtual void ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK)=0;

    /**
     * This function responsible for the  scalar variable projection calculation
     * @param elmtinfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param ctan 1x2 vector, where ctan[0] is responsible for the non-time-derivative part in the K matrix, while ctan[1] represents the coeffecient for the time derivatives in the K matrix
     * @param soln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param Mate the materials of current step
     * @param MateOld the materials of previous step
     * @param gpProj the projected scalar variable in current element
     */
    virtual void ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj)=0;

};
