//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.07.24
//+++ Purpose: the equation system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>

#include "MathUtils/Vector.h"
#include "MathUtils/SparseMatrix.h"
#include "DofHandler/BulkDofHandler.h"

/**
 * This class defines the sparse vector for the right hand side residual, and the spare matrix
 * for the system's jacobian matrix.
 */
class EquationSystem{
public:
    /**
     * constructor
     */
    EquationSystem();
    /**
     * initialize and allocate memory for residual and K matrix
     * @param t_dofHandler the dof handler class
     */
    void init(const DofHandler &t_dofHandler);
    /**
     * get the dofs number
     */
    inline int getDofsNum()const{return m_dofs;}

    /**
     * create the sparsity pattern for sparse matrix, after this call, 
     * the matrix structure can not be modified!!!
     * @param t_dofHandler the dof handler class
     */
    void createSparsityPattern(const DofHandler &t_dofHandler);
    /**
     * release the allocated memory
     */
    void releaseMemory();

public:
    Vector m_rhs;/**< vector for system residual */
    SparseMatrix m_amatrix;/**< sparse matrix for system K matrix */

private:
    int m_dofs;/**< the dofs or dimension of K matrix */
    bool m_allocated;/**< boolean flag for the memory allocation status */

};