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

#include "nlohmann/json.hpp"

#include "Utils/MessagePrinter.h"
#include "Utils/JsonUtils.h"

#include "MathUtils/Vector3d.h"
#include "MathUtils/VectorXd.h"
#include "MathUtils/Vector.h"
#include "MathUtils/MatrixXd.h"

#include "ElmtSystem/LocalElmtData.h"
#include "FESystem/FECalcType.h"

using std::string;
using std::vector;
using std::map;

/**
 * This class defines the abstract class for the integrated boundary condition
 */
class IntegrateBCBase{
public:
    /**
     * execute the boundary integration calculation for different integrated boundary conditions.
     * \f$\int_{\partial\Omega}f(u,\nabla u,\vec{n})dS\f$ is the final output
     * @param CalcType the calculation type for either residual or jacbobian
     * @param BCValue the boundary value defined in the input file
     * @param Params the json content for boundary condition parameters
     * @param ElmtInfo the basic information for current element
     * @param ElmtSoln the solution of current element
     * @param Normal the normal vector of current gauss point in the current element
     * @param Shp the local shape function values
     * @param Ctan the time derivative related coefficient array
     * @param LocalK the local K matrix
     * @param LocalR the locak R matrix
     */
    virtual void computeBCValue(const FECalcType &CalcType,
                                const double &BCValue,
                                const nlohmann::json Params,
                                const LocalElmtInfo &ElmtInfo,
                                const LocalElmtSolution &ElmtSoln,
                                const Vector3d &Normal,
                                const LocalShapeFun &Shp,
                                const double (&Ctan)[3],
                                MatrixXd &LocalK,
                                VectorXd &LocalR)=0;
    
protected:
    /**
     * calculate the residual of current boundary element
     * \f$\int_{\partial\Omega}R_{u}dS\f$ is the final output
     * @param BCValue the boundary value defined in the input file
     * @param Params the boundary condition related parameters (json content)
     * @param ElmtInfo the basic information for current element
     * @param ElmtSoln the solution of current element
     * @param Normal the normal vector of current gauss point in the current element
     * @param Shp the local shape function values
     * @param LocalR the locak R matrix
     */
    virtual void computeResidual(const double &BCValue,
                                 const nlohmann::json Params,
                                 const LocalElmtInfo &ElmtInfo,
                                 const LocalElmtSolution &ElmtSoln,
                                 const Vector3d &Normal,
                                 const LocalShapeFun &Shp,
                                 VectorXd &LocalR)=0; 
    
    /**
     * calculate the jacbobian contribution of current boundary element 
     * \f$\int_{\partial\Omega}\frac{\partial R{u}}{\partial u}dS\f$ is the final output
     * @param BCValue the boundary value defined in the input file
     * @param Params the boundary condition related parameters (json content)
     * @param ElmtInfo the basic information for current element
     * @param ElmtSoln the solution of current element
     * @param Normal the normal vector of current gauss point in the current element
     * @param Shp the local shape function values
     * @param LocalK the local K matrix
     */
    virtual void computeJacobian(const double &BCValue,
                                 const nlohmann::json &Params,
                                 const LocalElmtInfo &ElmtInfo,
                                 const LocalElmtSolution &ElmtSoln,
                                 const Vector3d &Normal,
                                 const LocalShapeFun &Shp,
                                 const double (&Ctan)[3],
                                 MatrixXd &LocalK)=0;

};