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
//+++ Date   : 2021.09.04
//+++ Purpose: define the abstract class for boundary condition with
//+++          integration on specific boundary, i.e., neumann bc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>

#include "Utils/MessagePrinter.h"

#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"
#include "Utils/MathFuns.h"

#include "ElmtSystem/LocalElmtData.h"
#include "FESystem/FECalcType.h"

using namespace std;

/**
 * This class defines the abstract class for integrated boundary condition
 */
class IntegrateBCBase{
public:
    /**
     * execute the boundary integration calculation for different integrated boundary conditions.
     * \f$\int_{\partial\Omega}f(u,\nabla u,\vec{n})dS\f$ is the final output
     * @param calctype the calculation type for either residual or jacbobian
     * @param bcvalue the boundary value defined in the input file
     * @param params the boundary condition related parameters
     * @param elmtinfo the basic information for current element
     * @param soln the solution of current element
     * @param normal the normal vector of current gauss point in the current element
     * @param shp the local shape function values
     * @param localK the local K matrix
     * @param localR the locak R matrix
     */
    virtual void ComputeBCValue(const FECalcType &calctype,const double &bcvalue,const vector<double> &params,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &soln,const Vector3d &normal,const LocalShapeFun &shp,MatrixXd &localK,VectorXd &localR)=0;
    
    /**
     * calculate the residual of current boundary element
     * \f$\int_{\partial\Omega}R_{u}dS\f$ is the final output
     * @param bcvalue the boundary value defined in the input file
     * @param params the boundary condition related parameters
     * @param elmtinfo the basic information for current element
     * @param soln the solution of current element
     * @param normal the normal vector of current gauss point in the current element
     * @param shp the local shape function values
     * @param localR the locak R matrix
     */
    virtual void ComputeResidual(const double &bcvalue,const vector<double> &params,const LocalElmtInfo &bcelmtinfo,const LocalElmtSolution &bcelmtsoln,const Vector3d &normal,const LocalShapeFun &shp,VectorXd &localR)=0; 
    
    /**
     * calculate the jacbobian contribution of current boundary element 
     * \f$\int_{\partial\Omega}\frac{\partial R{u}}{\partial u}dS\f$ is the final output
     * @param bcvalue the boundary value defined in the input file
     * @param params the boundary condition related parameters
     * @param elmtinfo the basic information for current element
     * @param soln the solution of current element
     * @param normal the normal vector of current gauss point in the current element
     * @param shp the local shape function values
     * @param localK the local K matrix
     */
    virtual void ComputeJacobian(const double &bcvalue,const vector<double> &params,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &soln,const Vector3d &normal,const LocalShapeFun &shp,MatrixXd &localK)=0;

};


