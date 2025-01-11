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
//+++ Date   : 2025.01.11
//+++ Purpose: This base class defines the base class for the linear
//+++          equation solver, which calculate Ax=b
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MathUtils/SparseMatrix.h"
#include "MathUtils/Vector.h"


class LinearSolverBase {
public:
    /**
     * Solve the linear equation Ax=b, where A is the sparse matrix, b is the right hand side vector, and x is the solution
     * @param A the input sparse matrix
     * @param b the rand hand side vector
     * @param x the solution vector
     * @return if success then return true, otherwise return false
     */
    bool solve(SparseMatrix &A,Vector &b,Vector &x)=0;

    /**
     * Get the iterations of the current linear equation solver
     * @return return the iteration number of the current linear equation solver
     */
    int getIterationNumber()const=0;
};