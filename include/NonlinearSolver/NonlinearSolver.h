//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_NONLINEARSOLVER_H
#define ASFEM_NONLINEARSOLVER_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

//************************************
//*** For AsFem's own header file
//************************************
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCs/BCSystem.h"
#include "ICs/ICSystem.h"

#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"

#include "Solution/Solution.h"
#include "EquationSystem/EquationSystem.h"
#include "FE/FE.h"
#include "FESystem/FESystem.h"


#include "NonlinearSolverBlock.h"
#include "NonlinearSolverType.h"
#include "FEProblem/FeCtrlInfo.h"

using namespace std;

typedef struct{
    Mesh _mesh;
    DofHandler _dofHandler;
    BCSystem _bcSystem;
    ICSystem _icSystem;
    ElmtSystem _elmtSystem;
    MateSystem _mateSystem;
    Solution _solution;
    EquationSystem _equationSystem;
    FE _fe;
    FESystem _feSystem;
    FeCtrlInfo _fectrlinfo;
} AppCtx;

typedef struct{
    PetscReal rnorm,rnorm0;
    PetscReal dunorm,dunorm0;
    PetscReal enorm,enorm0;
    PetscInt iters;
    bool IsDepDebug;
} MonitorCtx;


extern PetscErrorCode Monitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx);
extern PetscErrorCode FormJacobian(SNES snes,Vec U,Mat A,Mat B,void *ctx);
extern PetscErrorCode FormResidual(SNES snes,Vec U,Vec RHS,void *ctx);

extern PetscErrorCode MyConvergent(SNES snes,PetscInt iters,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason, void *cctx);

class NonlinearSolver{
public:
    NonlinearSolver();
    void Init(NonlinearSolverBlock &nonlinearsolverblock);

    inline PetscInt GetCurrentIters()const{return _Iters;}

    bool Solve(Mesh &mesh,DofHandler &dofHandler,
               ElmtSystem &elmtSystem,MateSystem &mateSystem,
               BCSystem &bcSystem,ICSystem &icSystem,
               Solution &solution,EquationSystem &equationSystem,
               FE &fe,FESystem &feSystem);
    bool SSolve(Mesh &mesh,DofHandler &dofHandler,
               ElmtSystem &elmtSystem,MateSystem &mateSystem,
               BCSystem &bcSystem,ICSystem &icSystem,
               Solution &solution,EquationSystem &equationSystem,
               FE &fe,FESystem &feSystem,
               FeCtrlInfo &fectrlinfo);

    inline double GetRnorm()const{return _Rnorm;}
    inline double GetdUnorm()const{return _dUnorm;}
    inline double GetEnorm()const{return _Enorm;}

private:
    bool NewtonRaphson(Mesh &mesh,DofHandler &dofHandler,
               ElmtSystem &elmtSystem,MateSystem &mateSystem,
               BCSystem &bcSystem,ICSystem &icSystem,
               Solution &solution,EquationSystem &equationSystem,
               FE &fe,FESystem &feSystem);
    

private:
    KSP _ksp;
    PC _pc;
    SNES _snes;
    SNESLineSearch _linesearch;
    SNESConvergedReason _snesreason;
    AppCtx _appctx;
    MonitorCtx _monctx;
    KSPConvergedReason _reason;

    bool LinearSolve(Mat &A,Vec &x,Vec &F);
    bool CheckConvergence();

    void PrintIterationInfo()const;
    void PrintIterationDetailsInfo()const;

private:
    //*********************************************
    //*** For nonlinear solver information
    //*********************************************
    PetscReal _Rnorm0,_Rnorm;
    PetscReal _Enorm0,_Enorm;
    PetscReal _dUnorm,_dUnorm0;
    PetscReal _RAbsTol,_RRelTol;
    PetscReal _EAbsTol,_ERelTol;
    PetscReal _STol;
    PetscInt _MaxIters,_Iters;
    bool _IsConvergent;
    NonlinearSolverType _SolverType;
};

#endif // ASFEM_NONLINEARSOLVER_H