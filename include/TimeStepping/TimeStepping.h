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
//+++ Date   : 2020.12.29
//+++ Purpose: Define time stepping system for transient analysis
//+++          in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>

#include "Utils/MessagePrinter.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"

#include "BCSystem/BCSystem.h"
#include "ICSystem/ICSystem.h"

#include "MateSystem/MateSystem.h"
#include "ElmtSystem/ElmtSystem.h"

#include "FE/FE.h"
#include "FESystem/FESystem.h"

#include "SolutionSystem/SolutionSystem.h"
#include "EquationSystem/EquationSystem.h"

#include "OutputSystem/OutputSystem.h"

#include "NonlinearSolver/NonlinearSolverBlock.h"

#include "TimeStepping/TimeSteppingBlock.h"
#include "TimeStepping/TimeSteppingType.h"

#include "FEProblem/FEControlInfo.h"

using namespace std;

//************************************************************************
//*** for some user-defined contex
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
    OutputSystem _outputSystem;
    FEControlInfo _fectrlinfo;
    //********************************
    PetscReal rnorm,rnorm0;
    PetscReal dunorm,dunorm0;
    PetscReal enorm,enorm0;
    PetscInt iters;
    bool IsDebug;
    bool IsDepDebug;
    double time;
    double dt;
    int step;
    //**************************
    bool IsAdaptive;
    int OptiIters;
    double GrowthFactor;
    double CutbackFactor;
    double DtMin;
    double DtMax;
} TSAppCtx;

//************************************************************************
//*** subroutines for the calculation of Residual and Jacobian
extern PetscErrorCode ComputeIResidual(TS ts,PetscReal t,Vec U,Vec V,Vec RHS,void *ctx);
extern PetscErrorCode ComputeIJacobian(TS ts,PetscReal t,Vec U,Vec V,PetscReal s,Mat A,Mat B,void *ctx);
extern PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal time,Vec U,void *ctx);
extern PetscErrorCode MySNESMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx);
//************************************************************************

class TimeStepping{
public:
    TimeStepping();

    void Init();

    void SetOpitonsFromTimeSteppingBlock(TimeSteppingBlock &timeSteppingBlock);

    void SetOptionsFromNonlinearSolverBlock(NonlinearSolverBlock &nonlinearsolverblock);

    bool Solve(Mesh &mesh,DofHandler &dofHandler,
            ElmtSystem &elmtSystem,MateSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            FE &fe,FESystem &feSystem,
            OutputSystem &outputSystem,
            FEControlInfo &fectrlinfo);

    void ReleaseMem();

    void PrintTimeSteppingInfo()const;
private:
    //*****************************************************************
    //*** basic variables for time stepping
    //*****************************************************************
    double _Dt=1.0e-5;
    double _FinalT=1.0e-3;
    bool _Adaptive=false;
    long int _TotalSteps=-1;
    long int _CurrentStep=-1;
    TimeSteppingType _TimeSteppingType=TimeSteppingType::BACKWARDEULER;
    string _TimeSteppingTypeName="backward-euler";
    double _GrowthFactor=1.1,_CutBackFactor=0.85;
    int _OptIters;
    double _DtMin,_DtMax;

private:
    //*****************************************************************
    //*** for nonlinear solver options
    //*****************************************************************
    double _Rnorm0,_Rnorm;
    double _Enorm0,_Enorm;
    double _dUnorm,_dUnorm0;
    double _RAbsTol,_RRelTol;
    double _EAbsTol,_ERelTol;
    double _STol;
    int _MaxIters,_Iters;
    bool _IsConvergent;
    NonlinearSolverType _SolverType;
    string _PCTypeName;
    //*****************************************************************
    //*** for TS components from PETSc
    //*****************************************************************
    TS _ts;
    TSAdapt _tsadapt;
    SNES _snes;
    KSP _ksp;
    PC _pc;
    SNESConvergedReason _snesreason;
    TSAppCtx _appctx;
};