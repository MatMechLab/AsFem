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
//+++ Date   : 2022.08.12
//+++ Purpose: the wrapper for SNES solver (PETSc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "LinearSolver/LinearSolver.h"
#include "NonlinearSolver/NonlinearSolverBase.h"
#include "NonlinearSolver/NonlinearSolverBlock.h"

/**
 * This is the struct that will be used to pass down/up the args we need to call SNES
 */
typedef struct{
    FECell *_fecell;
    DofHandler *_dofHandler;
    BCSystem *_bcSystem;
    ElmtSystem *_elmtSystem;
    MateSystem *_mateSystem;
    SolutionSystem *_solutionSystem;
    EquationSystem *_equationSystem;
    FE *_fe;
    FESystem *_feSystem;
    FEControlInfo *_fectrlinfo;
} AppCtx;

/**
 * The structure for SNES monitor
 */
typedef struct{
    double rnorm,rnorm0;
    double dunorm,dunorm0;
    double enorm,enorm0;
    int iters;
    bool IsDepDebug;
} MonitorCtx;

/**
 * This is the extern function from PETSc, you must keep it and offer an implementation
 * @param snes the SNES solver class
 * @param iters current iteration
 * @param rnorm the residual norm
 * @param ctx the user-defined content
 */
extern PetscErrorCode myMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx);

/***
 * This is the PETSc convergence monitor for SNES
 * @param snes the SNES solver class
 * @param iters current iteration
 * @param xnorm the solution norm
 * @param snorm the delta u norm
 * @param fnorm the residual norm
 * @param reason the converged reason
 * @param ctx the user-defined content
 */
extern PetscErrorCode myConvergence(SNES snes,PetscInt iters,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason, void *ctx);

/**
 * the jacobian evalution function
 * @param snes the SNES solver class
 * @param U the solution vector
 * @param A the jacobian matrix
 * @param B the precondtion matrix
 * @param ctx the user-defined content
 */
extern PetscErrorCode computeJacobian(SNES snes,Vec U,Mat A,Mat B,void *ctx);

/**
 * the residual evalution function
 * @param snes the SNES solver class
 * @param U the solution vector
 * @param RHS the residual vector
 * @param ctx the user-defined content
 */
extern PetscErrorCode computeResidual(SNES snes,Vec U,Vec RHS,void *ctx);

/**
 * This class implement the SNES wrapper based on PETSc package for nonlinear problem
 */
class SNESSolver:public NonlinearSolverBase{
public:
    SNESSolver();

    /**
     * setup the SNES solver from [nlsolver] block
     * @param t_nlblock the nonlinear solver block
     */
    void setFromNonlinearSolverBlock(const NonlinearSolverBlock &t_nlblock);

    /**
     * initialize the SNES solver
     * @param lsolver the linear solver class
     */
    void initSolver(LinearSolver &lsolver);
    /**
     * release the allocated memory in SNES solver
     */
    void releaseMemory();

    /**
     * solve the nonlinear equation, if success then return true
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dof class
     * @param t_FE the fe class
     * @param t_ElmtSystem the element system class
     * @param t_MateSystem the material system class
     * @param t_FESystem the fe system class
     * @param t_BCSystem the boundary condition system
     * @param t_SolnSystem the solution system class
     * @param t_EqSystem the equation system class
     * @param t_LinearSolver the linear solver system
     * @param t_FECtrlInfo the fe control info
     */
    virtual bool solve(FECell &t_FECell,
                       DofHandler &t_DofHandler,
                       FE &t_FE,
                       ElmtSystem &t_ElmtSystem,
                       MateSystem &t_MateSystem,
                       FESystem &t_FESystem,
                       BCSystem &t_BCSystem,
                       SolutionSystem &t_SolnSystem,
                       EquationSystem &t_EqSystem,
                       LinearSolver &t_LinearSolver,
                       FEControlInfo &t_FECtrlInfo) override;

    /**
     * get the nonlinear solver name
     */
    inline string getNonLinearSolverName()const{return m_NLSolverName;}

    /**
     * get the iteration number
     */
    inline int getIterationNum()const{return m_Iterations;}

    /**
     * get the initial norm of residual
     */
    inline double getInitResidualNorm()const{return m_RNorm0;}
    /**
     * get the residual norm
     */
    inline double getResidualNorm()const{return m_RNorm;}

    /**
     * print out the SNES solver information
     */
    void printSolverInfo()const;

private:
    bool m_Initialized;/**< boolean flag for the status of initializing */
    int m_MaxIters;/**< the maximum iterations */
    int m_Iterations;/**< the iteration number */

    double m_STol;/**< the tolerance for delta U */

    double m_AbsTolR;/**< the absolute tolerance for residual */
    double m_RelTolR;/**< the relative tolerance for residual */

    double m_AbsTolDu;/**< the absolute tolerance for delta u */
    double m_RelTolDu;/**< the relative tolerance for delta u */

    double m_AbsTolE;/**< the absolute tolerance for energy */
    double m_RelTolE;/**< the relative tolerance for energy */

    double m_RNorm0;/**< the initial norm of residual */
    double m_RNorm;/**< the intermediate or final norm of reisudal */

private:
    string m_NLSolverName;/**< the nonlinear solver name in SNES */
    NonlinearSolverType m_NLSolverType;/**< the nonlinear solver type */

private:
    SNES m_SNES;/**< the SNES solver from PETSc */
    SNESLineSearch m_SNESLineSearch;/**< for the line search ctx */
    SNESConvergedReason m_SNESConvergeReason;/**< for the converged reason ctx */

    AppCtx m_AppCtx;
    MonitorCtx m_MonCtx;

};