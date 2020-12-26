//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
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

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCSystem/BCSystem.h"
#include "ICSystem/ICSystem.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "FE/FE.h"
#include "FESystem/FESystem.h"
#include "EquationSystem/EquationSystem.h"
#include "NonlinearSolver/NonlinearSolverBlock.h"


//*******************************************************************

typedef struct{
    Mesh _mesh;
    DofHandler _dofHandler;
    BCSystem _bcSystem;
    ICSystem _icSystem;
    ElmtSystem _elmtSystem;
    MateSystem _mateSystem;
    SolutionSystem _solutionSystem;
    EquationSystem _equationSystem;
    FE _fe;
    FESystem _feSystem;
    // FeCtrlInfo _fectrlinfo;
} AppCtx;

typedef struct{
    PetscReal rnorm,rnorm0;
    PetscReal dunorm,dunorm0;
    PetscReal enorm,enorm0;
    PetscInt iters;
    bool IsDepDebug;
} MonitorCtx;

extern PetscErrorCode MyMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx);
extern PetscErrorCode MyConvergent(SNES snes,PetscInt iters,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason, void *cctx);

//************************************************************************
//*** the core part for our jacobian and residual calculation
//************************************************************************
extern PetscErrorCode FormJacobian(SNES snes,Vec U,Mat A,Mat B,void *ctx);
extern PetscErrorCode FormResidual(SNES snes,Vec U,Vec RHS,void *ctx);



class NonlinearSolver{
public:
    NonlinearSolver();
    void Init(NonlinearSolverBlock nonlinearsolverblock);
    bool Solve();

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
    string _SolverName;
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