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
//+++ Date   : 2020.07.12
//+++ Purpose: define the nonlinear solver class in AsFem
//+++          this class mainly call the SNES subroutines of PETSc
//+++          to solve the nonlinear equations
//+++          once PETSc update the API, we should also update the
//+++          related code !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

#include "Utils/MessagePrinter.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCSystem/BCSystem.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "FE/FE.h"
#include "FESystem/FESystem.h"
#include "EquationSystem/EquationSystem.h"
#include "NonlinearSolver/NonlinearSolverBlock.h"
#include "FEProblem/FEControlInfo.h"


//*******************************************************************
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
    PetscReal rnorm,rnorm0;
    PetscReal dunorm,dunorm0;
    PetscReal enorm,enorm0;
    PetscInt iters;
    bool IsDepDebug;
} MonitorCtx;

/**
 * This is the extern function from PETSc, you must keep it and offer an implementation
 */
extern PetscErrorCode MyMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx);

/***
 * This is the PETSc convergence monitor for SNES
 */
extern PetscErrorCode MyConvergent(SNES snes,PetscInt iters,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason, void *cctx);

/**
 * the core part for our jacobian and residual calculation
 */
extern PetscErrorCode ComputeJacobian(SNES snes,Vec U,Mat A,Mat B,void *ctx);
extern PetscErrorCode ComputeResidual(SNES snes,Vec U,Vec RHS,void *ctx);


/**
 * The nonlinear solver class, where we do the API wrapper for PETSc/SNES
 */
class NonlinearSolver{
public:
    /**
     * the constructor for nonlinear solver
     */
    NonlinearSolver();
    
    /**
     * Initial settings for nonlinear solver
     */
    void Init();
   
    /**
     * Get the final iterations of current solution
     */
    int GetFinalInterations()const{return _Iters;}  
    
    /**
     * Get the initial Residual norm of current solution
     */
    double GetInitialRNorm()const{return _Rnorm0;}

    /**
     * Get the final Residual norm of current solution
     */
    double GetFinalRNorm()const{return _Rnorm;}

    /**
     * Do the options/settings which are read from the [nonlinearsolver] block
     */
    void SetOptionsFromNonlinearSolverBlock(NonlinearSolverBlock &nonlinearsolverblock);
    
    /**
     * Solve the nonlinear equation via SNES, if the system is converged then return true else the  false is returned.
     * @param mesh the mesh class
     * @param dofHandler the dof manager class
     * @param elmtSystem the elmtSystem class
     * @param mateSystem the material system class
     * @param bcSystem the boundary condition class
     * @param icSystem the initial condition class
     * @param solutionSystem the solution class
     * @param equationSystem the equation system class
     * @param fe the FE space class
     * @param feSystem the feSystem for FE calculation
     * @param fectrlinfo the fe control structure
     */
    bool Solve(Mesh &mesh,DofHandler &dofHandler,
            ElmtSystem &elmtSystem,MateSystem &mateSystem,
            BCSystem &bcSystem,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            FE &fe,FESystem &feSystem,
            FEControlInfo &fectrlinfo);

    /**
     * Release the allocated memory
     */
    void ReleaseMem();

    /**
     * Print the basic information of nonlinearsolver
     */
    void PrintInfo()const;
private:
    //*********************************************
    //*** For nonlinear solver information
    //*********************************************
    double _Rnorm0,_Rnorm;
    double _Enorm0,_Enorm;
    double _dUnorm,_dUnorm0;
    double _RAbsTol,_RRelTol;
    double _EAbsTol,_ERelTol;
    double _STol;
    int _MaxIters,_Iters;
    bool _IsConvergent;
    NonlinearSolverType _SolverType;
    string _LinearSolverName,_SolverTypeName;
    string _PCTypeName;

    //*********************************************
    //*** For nonlinear solver's related components
    //*********************************************
    KSP _ksp;
    PC  _pc;
    SNES _snes;
    SNESLineSearch _sneslinesearch;
    SNESConvergedReason _snesreason;
    AppCtx _appctx;
    MonitorCtx _monctx;

};
