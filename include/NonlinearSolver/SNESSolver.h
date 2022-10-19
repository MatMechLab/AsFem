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
//+++ Date   : 2022.08.12
//+++ Purpose: the wrapper for SNES solver (PETSc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "NonlinearSolver/NonlinearSolverBase.h"
#include "NonlinearSolver/NonlinearSolverBlock.h"

/**
 * This is the struct that will be used to pass down/up the args we need to call SNES
 */
typedef struct{
    Mesh *_mesh;
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
     */
    void initSolver();
    /**
     * release the allocated memory in SNES solver
     */
    void releaseMemory();

    /**
     * solve the nonlinear equation, if success then return true
     * @param t_mesh the mesh class
     * @param t_dofhandler the dof class
     * @param t_fe the fe class
     * @param t_elmtsyste the element system class
     * @param t_matesystem the material system class
     * @param t_fesystem the fe system class
     * @param t_bcsystem the boundary condition system
     * @param t_solutionsystem the solution system class
     * @param t_equationsystem the equation system class
     * @param t_fectrlinfo the fe control info
     */
    virtual bool solve(Mesh &t_mesh,DofHandler &t_dofhandler,FE &t_fe,
                       ElmtSystem &t_elmtsyste,MateSystem &t_matesystem,
                       FESystem &t_fesystem,
                       BCSystem &t_bcsystem,
                       SolutionSystem &t_solutionsystem,
                       EquationSystem &t_equationsystem,
                       FEControlInfo &t_fectrlinfo) override;
    
    /**
     * get the linear solver name in current SNES solver
     */
    inline string getLinearSolverName()const{return m_linearsolvername;}

    /**
     * get the nonlinear solver name
     */
    inline string getNonLinearSolverName()const{return m_nlsolvername;}

    /**
     * get the iteration number
     */
    inline int getIterationNum()const{return m_iterations;}

    /**
     * get the initial norm of residual
     */
    inline double getInitResidualNorm()const{return m_rnorm0;}
    /**
     * get the residual norm
     */
    inline double getResidualNorm()const{return m_rnorm;}

    /**
     * print out the SNES solver information
     */
    void printSolverInfo()const;

private:
    bool m_initialized;/**< boolean flag for the status of initializing */
    int m_maxiters;/**< the maximum iterations */
    int m_iterations;/**< the iteration number */

    double m_s_tol;/**< the tolerance for delta U */

    double m_abstol_r;/**< the absolute tolerance for residual */
    double m_reltol_r;/**< the relative tolerance for residual */

    double m_abstol_du;/**< the absolute tolerance for delta u */
    double m_reltol_du;/**< the relative tolerance for delta u */

    double m_abstol_e;/**< the absolute tolerance for energy */
    double m_reltol_e;/**< the relative tolerance for energy */

    double m_rnorm0;/**< the initial norm of residual */
    double m_rnorm;/**< the intermediate or final norm of reisudal */

private:
    string m_linearsolvername;/**< the string name of the linear solver in SNES*/
    string m_nlsolvername;/**< the nonlinear solver name in SNES */
    string m_pcname;/**< the preconditioner name of current SNES solver */
    NonlinearSolverType m_nlsolvertype;/**< the nonlinear solver type */

private:
    SNES m_snes;/**< the SNES solver from PETSc */
    KSP  m_ksp;/**< the KSP linear solver from PETSc */
    PC   m_pc;/**< the preconditioner from PETSc */
    SNESLineSearch m_sneslinesearch;/**< for the line search ctx */
    SNESConvergedReason m_snesconvergereason;/**< for the converged reason ctx */

    AppCtx m_appctx;
    MonitorCtx m_monctx;

};