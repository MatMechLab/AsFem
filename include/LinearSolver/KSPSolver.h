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
//+++ Purpose: This class offers the API for the KSP linear solver
//+++          from PETSc, one should have PETSc installed
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "LinearSolver/LinearSolverBase.h"

#include "petsc.h"

class KSPSolver:public LinearSolverBase {
public:
    /**
     * constructor
     */
    KSPSolver();

    /**
     *
     * @param A the SparseMatrix
     * @param b the right hand side vector
     * @param x the solution vector
     * @return if solve success return true, otherwise return false
     */
    virtual bool solve(SparseMatrix &A,Vector &b,Vector &x) override;

    /**
     * init the linear solver
     */
    void init();

    /**
     * release the allocated memory
     */
    void releaseMemory();

    /**
     * set the default parameters for KSP solver
     */
    void setDefaultParams();

    /**
     * set the GMRES restart number of ksp solver
     * @param restartNumber the restart number
     */
    void setGMRESRestartNumber(const int &restartNumber) {
        m_GMRESRestartNumber = restartNumber;
    }

    /**
     * set the maximum iterations of the linar solver
     * @param maxIterations the maximum iteration for the linear solver
     */
    void setMaxIterations(const int &maxIterations) {
        m_MaxIterations=maxIterations;
    }
    /**
     * set the tolerance of ksp solver
     * @param tolerance the tolerance of ksp solver
     */
    void setSolverTolerance(const double &tolerance) {
        m_Tolerance = tolerance;
    }

    /**
     * set the solver name of ksp
     * @param typeName string name for the ksp solver
     */
    void setKSPSolverTypeName(const std::string &typeName) {
        m_KSPSolverTypeName = typeName;
    }

    /**
     * set the preconditioner name in ksp
     * @param typeName the string name for the preconditioner in KSP
     */
    void setKSPPCTypeName(const std::string &typeName) {
        m_PCTypeName = typeName;
    }

    /**
     * get the KSP class copy
     * @return the copy of KSP class
     */
    KSP getKSPCopy()const {
        return m_KSP;
    }
    /**
     * get the reference of the KSP class
     * @return the reference of the KSP class
     */
    KSP& getKSPRef() {
        return m_KSP;
    }

    /**
     * get the copy of PC class
     * @return PC copy
     */
    PC getPCCopy()const {
        return m_PC;
    }
    /**
     * get the reference of PC class
     * @return PC reference
     */
    PC& getPCRef() {
        return m_PC;
    }

    /**
     * get the iteration number
     * @return iteration number
     */
    int getIterationNumber()const {
        return m_Iterations;
    }

    /**
     * get the restart number of GMRES solver
     * @return GMRES restart number
     */
    int getGMRESRestartNumber()const {
        return m_GMRESRestartNumber;
    }

    /**
     * get the solver name of ksp class
     * @return the string name of ksp solver
     */
    string getSolverTypeName()const {
        return m_KSPSolverTypeName;
    }

    /**
     * get the preconditoner name of ksp's pc
     * @return the string name of pc type
     */
    string getPCTypeName()const {
        return m_PCTypeName;
    }

    /**
     * print the summary information of KSP solver
     */
    void printKSPSolverInfo()const;


protected:
    int m_Iterations;/**< the iterations of the linear solver */
    int m_MaxIterations;/**< the maximum iterations of the linear solver */
    int m_GMRESRestartNumber;/**< the restart number of GMRES */
    double m_Tolerance;/**< the error tolerance of the linear solver */
    bool m_IsAllocated;/**< the status of the memory allocation */

    KSP m_KSP;/**< the KSP solver class */
    PC m_PC;/**< the Preconditoner class */

    string m_KSPSolverTypeName;/**< the string name of the KSP Solver type */
    string m_PCTypeName;/**< the string name of the precondition type */
};
