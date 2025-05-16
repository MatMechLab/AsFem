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
//+++ Date   : 2022.08.26
//+++ Purpose: implement the dirichlet boundary condition for 2d
//+++          poisson benchmark test
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "BCSystem/DirichletBCBase.h"

/**
 * This class implement the dirichelt boundary condition calculation for 2d poisson benchmark test
 */
class Poisson2DBenchmarkBC:public DirichletBCBase{
public:
    Poisson2DBenchmarkBC();
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
                                Vector &RHS) override;

private:
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
                          VectorXd &LocalU) override; 

};
