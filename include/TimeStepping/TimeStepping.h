//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_TIMESTEPPER_H
#define ASFEM_TIMESTEPPER_H

#include <iostream>
#include <iomanip>
#include <fstream>

#include "petsc.h"

//***************************************
//*** For AsFem's own header file
//***************************************
#include "TimeSteppingType.h"
#include "TimeSteppingBlock.h"
#include "NonlinearSolver/NonlinearSolverBlock.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCs/BCSystem.h"
#include "ICs/ICSystem.h"

#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"

#include "EquationSystem/EquationSystem.h"
#include "Solution/Solution.h"

#include "FE/FE.h"
#include "FESystem/FESystem.h"

#include "OutputSystem/OutputSystem.h"

#include "NonlinearSolver/NonlinearSolver.h"

#include "FEProblem/FeCtrlInfo.h"

using namespace std;


//***************************************
//**** basic components for PETSc
//***************************************
// extern PetscErrorCode Residual(TS ts,PetscReal t,Vec U,Vec Udot,Vec RHS,void *ctx);
// extern PetscErrorCode Jacobian(TS ts,PetscReal t,Vec U,Vec Udot,PetscReal s,Mat A,Mat B,void *ctx);

// extern PetscErrorCode MySNESMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx);
// extern PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal time,Vec U,void *ctx);

// typedef struct{
//     Mesh _mesh;
//     DofHandler _dofHandler;
//     BCSystem _bcSystem;
//     ICSystem _icSystem;
//     ElmtSystem _elmtSystem;
//     MateSystem _mateSystem;
//     Solution _solution;
//     EquationSystem _equationSystem;
//     FE _fe;
//     FESystem _feSystem;
//     OutputSystem _outputSystem;
//     //******************************
//     PetscReal rnorm,rnorm0;
//     PetscReal dunorm,dunorm0;
//     PetscReal enorm,enorm0;
//     PetscReal time,dt;
//     PetscInt iters;
//     PetscInt step;
//     PetscInt interval;
//     bool IsDebug;
//     bool IsDepDebug;
//     bool IsProjection;
// } TSAppCtx;




class TimeStepping{
public:
    TimeStepping();
    void Init(TimeSteppingBlock &timesteppingblock);
    // void InitSolver(NonlinearSolverBlock &nonlinearsolverblock);

    void Stepping(/*Mesh &mesh,DofHandler &dofHandler,
                BCSystem &bcSystem,ICSystem &icSystem,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                EquationSystem &equationSystem,Solution &solution,
                FE &fe,FESystem &feSystem,
                OutputSystem &outputSystem,
                FeCtrlInfo &fectrl*/);
    
    void SteppingNew(Mesh &mesh,DofHandler &dofHandler,
                BCSystem &bcSystem,ICSystem &icSystem,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                EquationSystem &equationSystem,Solution &solution,
                FE &fe,FESystem &feSystem,
                OutputSystem &outputSystem,
                NonlinearSolver &nonlinearsolver,
                FeCtrlInfo &fectrl);


    void PrintInfo() const;

private:
    void BackwardEuler(Mesh &mesh,DofHandler &dofHandler,
                BCSystem &bcSystem,ICSystem &icSystem,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                EquationSystem &equationSystem,Solution &solution,
                FE &fe,FESystem &feSystem,
                OutputSystem &outputSystem,
                NonlinearSolver &nonlinearsolver,
                FeCtrlInfo &fectrl);
        
    void CrankNicolson(Mesh &mesh,DofHandler &dofHandler,
                BCSystem &bcSystem,ICSystem &icSystem,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                EquationSystem &equationSystem,Solution &solution,
                FE &fe,FESystem &feSystem,
                OutputSystem &outputSystem,
                NonlinearSolver &nonlinearsolver,
                FeCtrlInfo &fectrl);


private:
    TS _ts;
    TSAdapt _adapt;
    SNES _snes;
    KSP _ksp;
    PC _pc;
    SNESLineSearch _linesearch;
    PetscReal _FinalTime;
    PetscReal _dtmax,_dtmin,_dt0;
    PetscInt _interval;
    PetscInt _nOpts=5;
    TimeSteppingType _TimeSteppingType;
    bool _IsAdaptive=false;
    PetscReal _CutFactor,_GrowthFactor;

    // TSAppCtx _appctx;

};

#endif // ASFEM_TIMESTEPPER_H