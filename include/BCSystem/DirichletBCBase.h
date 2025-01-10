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
//+++ Date   : 2021.10.06
//+++ Purpose: define the abstract class for dirichlet type boundary
//+++          condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>

#include "Utils/MessagePrinter.h"
#include "Utils/JsonUtils.h"

#include "MathUtils/Vector3d.h"
#include "MathUtils/VectorXd.h"
#include "MathUtils/Vector.h"
#include "MathUtils/MatrixXd.h"
#include "MathUtils/SparseMatrix.h"

#include "ElmtSystem/LocalElmtData.h"
#include "FESystem/FECalcType.h"

#include "nlohmann/json.hpp"

using std::string;
using std::vector;
using std::map;

/**
 * This class defines the abstract class for dirichlet type boundary condition
 */
class DirichletBCBase{
public:
    /**
     * execute the boundary condition value for different (dirichlet type) boundary conditions.
     * \f$u=u_{g}\f$ is the final output
     * @param CalcType the calculation type for either residual or jacbobian
     * @param Penalty the penalty for dirichlet boundary conditions
     * @param BCValue the boundary value defined in the input file
     * @param Params the boundary condition related parameters(json file content)
     * @param ElmtInfo the basic information for current element
     * @param ElmtSoln the element solution of current element
     * @param DofIDs the global id of the applied dof(start from 1, the global one)
     * @param U the system solution
     * @param K the system sparse matrix
     * @param RHS the system residual vector
     */
    virtual void computeBCValue(const FECalcType &CalcType,
                                const double &Penalty,
                                const double &BCValue,
                                const nlohmann::json &Params,
                                const LocalElmtInfo &ElmtInfo,
                                const LocalElmtSolution &ElmtSoln,
                                const vector<int> &DofIDs,
                                Vector &U,
                                SparseMatrix &K,
                                Vector &RHS)=0;

    /**
     * calculate the 'displacement' value of current dofs
     * @param BCValue the boundary value defined in the input file
     * @param Params the boundary condition related parameters(json file content)
     * @param DofIDs the global id of the applied dof(start from 1, the global one)
     * @param ElmtInfo the basic information for current element
     * @param ElmtSoln the solution of current element
     * @param LocalU the solution vector of current node
     */
    virtual void computeU(const double &BCValue,
                          const nlohmann::json &Params,
                          const vector<int> &DofIDs,
                          const LocalElmtInfo &ElmtInfo,
                          const LocalElmtSolution &ElmtSoln,
                          VectorXd &LocalU)=0;

protected:
    VectorXd m_LocalU;/**< local solution vector */

};

