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
//+++ Date   : 2020.12.30
//+++ Purpose: Create the TS solver for transient analysis in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TimeStepping/TimeStepping.h"

void TimeStepping::Init(){
    TSCreate(PETSC_COMM_WORLD,&_ts);
    TSSetProblemType(_ts,TS_NONLINEAR);// in AsFem, we only consider the implicit problem
//    TSSetEquationType(_ts,TS_EQ_IMPLICIT);
    TSSetEquationType(_ts,TS_EQ_DAE_IMPLICIT_INDEX3);
//    TSSetEquationType(_ts,TS_EQ_ODE_IMPLICIT);
    //***************************************************
    //*** for the different time stepping methods
    //***************************************************
    if(_TimeSteppingType==TimeSteppingType::BACKWARDEULER){
        TSSetType(_ts,TSBEULER);
    }
    else if(_TimeSteppingType==TimeSteppingType::CRANCKNICLSON){
        TSSetType(_ts,TSCN);
    }
    else if(_TimeSteppingType==TimeSteppingType::ALPHA){
        TSSetType(_ts,TSALPHA);
    }
    else if(_TimeSteppingType==TimeSteppingType::GL){
        TSSetType(_ts,TSGLLE);
    }
    else if(_TimeSteppingType==TimeSteppingType::ROSW){
        TSSetType(_ts,TSROSW);
    }
    else if(_TimeSteppingType==TimeSteppingType::BDF2){
        TSSetType(_ts,TSBDF);
        TSBDFSetOrder(_ts,2);
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported time stepping method");
        MessagePrinter::AsFem_Exit();
    }

    TSSetTimeStep(_ts,_Dt);// we set the initial delta t
    TSSetMaxTime(_ts,_FinalT);
    TSSetExactFinalTime(_ts,TS_EXACTFINALTIME_INTERPOLATE);

    TSGetAdapt(_ts,&_tsadapt);
    TSAdaptSetType(_tsadapt,TSADAPTNONE);// disable adaptivity, we will use our own adapter
    TSSetTolerances(_ts,_RAbsTol,NULL,_RRelTol,NULL);
    cout<<"abs tol="<<_RAbsTol<<", rel tol="<<_RRelTol<<endl;
    //***************************************************
    //*** for the nonlinear solver settings
    //***************************************************
    TSGetSNES(_ts,&_snes);
    //**************************************************
    //*** init KSP
    //**************************************************
    SNESGetKSP(_snes,&_ksp);
    KSPGMRESSetRestart(_ksp,1400);
    KSPGetPC(_ksp,&_pc);
    PCFactorSetMatSolverType(_pc,MATSOLVERPETSC);

    if(_LinearSolverName=="mumps"){
        PCSetType(_pc,PCLU);
        KSPSetType(_ksp,KSPPREONLY);
        PCFactorSetMatSolverType(_pc,MATSOLVERMUMPS);
    }
    else if(_LinearSolverName=="superlu"){
        PCSetType(_pc,PCLU);
        KSPSetType(_ksp,KSPPREONLY);
        PCFactorSetMatSolverType(_pc,MATSOLVERSUPERLU_DIST);
    }
    else{
        PCSetType(_pc,PCLU);
    }

    PCFactorSetReuseOrdering(_pc,PETSC_TRUE);

    //**************************************************
    //*** allow user setting ksp from command line
    //**************************************************
    KSPSetFromOptions(_ksp);

    //**************************************************
    //*** some basic settings for SNES
    //**************************************************
    SNESSetTolerances(_snes,_RAbsTol,_RRelTol,_STol,_MaxIters,-1);
    SNESSetDivergenceTolerance(_snes,-1);
    TSSetMaxSNESFailures(_ts,-1);
    //**************************************************
    //*** for different type of nonlinear methods
    //**************************************************
    SNESSetType(_snes,SNESNEWTONLS);// our default method
    if(_SolverType==NonlinearSolverType::NEWTON||_SolverType==NonlinearSolverType::NEWTONLS){
        SNESSetType(_snes,SNESNEWTONLS);
    }
    else if(_SolverType==NonlinearSolverType::NEWTONTR){
        SNESSetType(_snes,SNESNEWTONTR);
    }
    else if(_SolverType==NonlinearSolverType::BFGS){
        SNESSetType(_snes,SNESQN);
    }
    else if(_SolverType==NonlinearSolverType::BROYDEN){
        SNESSetType(_snes,SNESQN);
        SNESQNSetType(_snes,SNES_QN_BROYDEN);
    }
    else if(_SolverType==NonlinearSolverType::BADBROYDEN){
        SNESSetType(_snes,SNESQN);
        SNESQNSetType(_snes,SNES_QN_BADBROYDEN);
    }
    else if(_SolverType==NonlinearSolverType::NEWTONCG){
        SNESSetType(_snes,SNESNCG);
    }
    else if(_SolverType==NonlinearSolverType::NEWTONGMRES){
        SNESSetType(_snes,SNESNGMRES);
    }

}
