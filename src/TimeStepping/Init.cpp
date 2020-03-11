//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "TimeStepping/TimeStepping.h"

void TimeStepping::Init(TimeSteppingBlock &timesteppingblock){
    _FinalTime=timesteppingblock._FinalTime;
    _dtmax=timesteppingblock._dtmax;
    _dtmin=timesteppingblock._dtmin;
    _dt0=timesteppingblock._dt0;
    _interval=timesteppingblock._interval;
    _TimeSteppingType=timesteppingblock._TimeSteppingMethod;
    _IsAdaptive=timesteppingblock._IsAdaptive;
    _nOpts=timesteppingblock._nOpts;
    _CutFactor=timesteppingblock._CutFactor;
    _GrowthFactor=timesteppingblock._GrowthFactor;
}
//****************************************************
//*** Now we init both the time stepping solver 
//*** and also the SNES solver
//****************************************************
// void TimeStepping::InitSolver(NonlinearSolverBlock &nonlinearsolverblock){
//     TSCreate(PETSC_COMM_WORLD,&_ts);
    
//     TSGetSNES(_ts,&_snes);
//     SNESSetForceIteration(_snes,PETSC_TRUE);
    
//     if(nonlinearsolverblock._SolverType==NonlinearSolverType::NewtonRaphson){
//         SNESSetType(_snes,SNESNEWTONLS);
//         SNESGetLineSearch(_snes,&_linesearch);
//         // SNESSetLineSearch(_snes,SNESLINESEARCHBASIC);
//         SNESLineSearchSetType(_linesearch,SNESLINESEARCHBASIC);
//     }
//     else if(nonlinearsolverblock._SolverType==NonlinearSolverType::SNESNewtonLs){
//         SNESSetType(_snes,SNESNEWTONLS);
//     }
//     else if(nonlinearsolverblock._SolverType==NonlinearSolverType::SNESNewtonTr){
//         SNESSetType(_snes,SNESNEWTONTR);
//     }
//     else if(nonlinearsolverblock._SolverType==NonlinearSolverType::SNESNewtonGMRES){
//         SNESSetType(_snes,SNESNGMRES);
//     }
//     else if(nonlinearsolverblock._SolverType==NonlinearSolverType::SNESNewtonCG){
//         SNESSetType(_snes,SNESNCG);
//     }
//     else if(nonlinearsolverblock._SolverType==NonlinearSolverType::SNESLBfgs){
//         SNESSetType(_snes,SNESQN);
//         SNESQNSetType(_snes,SNES_QN_LBFGS);
//     }
//     else if(nonlinearsolverblock._SolverType==NonlinearSolverType::SNESBroyden){
//         SNESSetType(_snes,SNESQN);
//         SNESQNSetType(_snes,SNES_QN_BROYDEN);
//     }
//     else if(nonlinearsolverblock._SolverType==NonlinearSolverType::SNESBadBroyden){
//         SNESSetType(_snes,SNESQN);
//         SNESQNSetType(_snes,SNES_QN_BADBROYDEN);
//     }
//     //**********************************
//     //*** for nonlinear solver tolerence settings
//     // SNESSetTolerances(_snes,
//     //                   nonlinearsolverblock._RAbsTol,
//     //                   nonlinearsolverblock._RRelTol,
//     //                   nonlinearsolverblock._STol,
//     //                   nonlinearsolverblock._MaxIters,
//     //                   -1);
//     SNESSetTolerances(_snes,
//                       nonlinearsolverblock._RAbsTol,
//                       nonlinearsolverblock._RRelTol,
//                       PETSC_DEFAULT,
//                       nonlinearsolverblock._MaxIters,
//                       -1);
//     SNESSetDivergenceTolerance(_snes,-1);
    
    
//     //**************************************
//     //*** Now we set the linear solver
//     //**************************************
//     SNESGetKSP(_snes,&_ksp);
//     KSPGetPC(_ksp,&_pc);
//     //KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,PetscReal dtol,PetscInt maxits)
//     KSPSetTolerances(_ksp,1.0e-10,1.0e-9,PETSC_DEFAULT,500000);
//     KSPGMRESSetRestart(_ksp,1200);
//     PCSetType(_pc,PCLU);

//     //*********************************************
//     //*** set the time stepping method
//     //*** and also the dt, dtmin, dtmax info
//     //*********************************************
//     TSSetProblemType(_ts,TS_NONLINEAR);
//     TSSetEquationType(_ts,TS_EQ_IMPLICIT);
//     if(_TimeSteppingType==TimeSteppingType::BackWardEuler){
//         TSSetType(_ts,TSBEULER);
//     }
//     else if(_TimeSteppingType==TimeSteppingType::CrankNicolson){
//         TSSetType(_ts,TSCN);
//     }
//     else if(_TimeSteppingType==TimeSteppingType::ThetaMethod){
//         TSSetType(_ts,TSTHETA);
//         TSThetaSetTheta(_ts,0.6);
//     }
//     else if(_TimeSteppingType==TimeSteppingType::AlphaMethod){
//         TSSetType(_ts,TSALPHA);
//         // TSAlphaSetParams()
//     }

//     // TSSetTime(_ts,_FinalTime);
//     TSSetMaxTime(_ts,_FinalTime);
//     TSSetTimeStep(_ts,_dt0);
//     TSSetMaxSteps(_ts,10000000);
//     if(_IsAdaptive){
//         TSGetAdapt(_ts,&_adapt);
//         TSAdaptSetStepLimits(_adapt,_dtmin,_dtmax);
//         TSAdaptSetType(_adapt,TSADAPTBASIC);
//         // TSAdaptSetType(_adapt,TSADAPTGLEE);
//         // TSAdaptSetType(_adapt,TSADAPTCFL);
//         TSAdaptSetClip(_adapt,0.8,1.1);
//         TSAdaptSetSafety(_adapt,0.9,0.85);
//         // TSAdaptSetTimeStepIncreaseDelay(_adapt,1);
//     }
//     else{
//         TSGetAdapt(_ts,&_adapt);
//         TSAdaptSetStepLimits(_adapt,_dtmin,_dt0);
//     }

    
//     TSSetExactFinalTime(_ts,TS_EXACTFINALTIME_INTERPOLATE);

//     TSSetMaxSNESFailures(_ts,-1);
//     // TSSetMaxStepRejections(_ts,6);
    

//     TSSetFromOptions(_ts);
//     SNESSetFromOptions(_snes);
//     KSPSetFromOptions(_ksp); 
//     PCSetFromOptions(_pc);

//     // TSSetUp(_ts);
//     // SNESSetUp(_snes);
    
    
// }