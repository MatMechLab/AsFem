//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
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

#include "NonlinearSolver/NonlinearSolver.h"

NonlinearSolver::NonlinearSolver(){
    _Rnorm0=1.0;_Rnorm=1.0;
    _Enorm0=1.0;_Enorm=1.0;
    _RAbsTol=4.0e-8;_RRelTol=1.0e-9;
    _EAbsTol=1.0e-19;_ERelTol=1.0e-20;
    _MaxIters=25;_Iters=0;
    _STol=1.0e-16;
    _SolverType=NonlinearSolverType::NEWTONLS;
    _SolverTypeName="newton with line search";
    _LinearSolverName="gmres";
    _PCTypeName="lu";
    _CheckJacobian=false;
}

void NonlinearSolver::SetOptionsFromNonlinearSolverBlock(NonlinearSolverBlock &nonlinearsolverblock){
    _SolverType=nonlinearsolverblock._SolverType;
    _MaxIters=nonlinearsolverblock._MaxIters;
    _RAbsTol=nonlinearsolverblock._RAbsTol;
    _RRelTol=nonlinearsolverblock._RRelTol;
    _STol=nonlinearsolverblock._STol;

    _SolverType=nonlinearsolverblock._SolverType;
    _SolverTypeName=nonlinearsolverblock._SolverTypeName;
    _LinearSolverName=nonlinearsolverblock._LinearSolverName;
    _PCTypeName=nonlinearsolverblock._PCTypeName;

    _CheckJacobian=nonlinearsolverblock._CheckJacobian;
}
void NonlinearSolver::Init(){
    //**************************************************
    //*** create our SNES solver
    //**************************************************
    SNESCreate(PETSC_COMM_WORLD,&_snes);

    //**************************************************
    //*** init KSP
    //**************************************************
    SNESGetKSP(_snes,&_ksp);
    KSPGMRESSetRestart(_ksp,1800);
    KSPGetPC(_ksp,&_pc);
    //PCFactorSetMatSolverType(_pc,MATSOLVERPETSC); // this line may be not so necessary!
    

    if(_LinearSolverName=="gmres"){
        KSPSetType(_ksp,KSPGMRES);
    }
    else if(_LinearSolverName=="fgmres"){
        KSPSetType(_ksp,KSPFGMRES);
    }
    else if(_LinearSolverName=="cg"){
        KSPSetType(_ksp,KSPCG);
    }
    else if(_LinearSolverName=="bicg"){
        KSPSetType(_ksp,KSPBICG);
    }
    else if(_LinearSolverName=="richardson"){
        KSPSetType(_ksp,KSPRICHARDSON);
    }
    else if(_LinearSolverName=="mumps"){
        PCSetType(_pc,PCLU);
        KSPSetType(_ksp,KSPPREONLY);
        PCFactorSetMatSolverType(_pc,MATSOLVERMUMPS);
    }
    else if(_LinearSolverName=="superlu"){
        PCSetType(_pc,PCLU);
        KSPSetType(_ksp,KSPPREONLY);
        PCFactorSetMatSolverType(_pc,MATSOLVERSUPERLU_DIST);
    }

    PCFactorSetReuseOrdering(_pc,PETSC_TRUE);

    //**************************************************
    //*** allow user setting ksp from command line
    //**************************************************
    KSPSetFromOptions(_ksp);
    PCSetFromOptions(_pc);
    //**************************************************
    //*** some basic settings for SNES
    //**************************************************
    SNESSetTolerances(_snes,_RAbsTol,_RRelTol,_STol,_MaxIters,-1);
    SNESSetDivergenceTolerance(_snes,-1);

    //**************************************************
    //*** for different type of nonlinear methods
    //**************************************************
    SNESSetType(_snes,SNESNEWTONLS);// our default method
    if(_SolverType==NonlinearSolverType::NEWTON || _SolverType==NonlinearSolverType::NEWTONLS){
        SNESSetType(_snes,SNESNEWTONLS);
        SNESGetLineSearch(_snes,&_sneslinesearch);
        SNESLineSearchSetType(_sneslinesearch,SNESLINESEARCHBT);
        SNESLineSearchSetOrder(_sneslinesearch,2);
    }
    else if(_SolverType==NonlinearSolverType::NEWTONSECANT){
        SNESSetType(_snes,SNESNEWTONLS);
        SNESGetLineSearch(_snes,&_sneslinesearch);
        SNESLineSearchSetType(_sneslinesearch,SNESLINESEARCHL2);
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

    SNESSetFromOptions(_snes);
}

//***************************************************
void NonlinearSolver::ReleaseMem(){
    SNESDestroy(&_snes);
}

//****************************************************
void NonlinearSolver::PrintInfo()const{
    MessagePrinter::PrintNormalTxt("Nonlinear solver information summary:");
    char buff[70];
    string str;

    if(_SolverType==NonlinearSolverType::NEWTON||_SolverType==NonlinearSolverType::NEWTONLS){
        str="  solver type= newton with line search";
    }
    else if(_SolverType==NonlinearSolverType::NEWTONTR){
        str="  solver type= newton trust region";
    }
    else if(_SolverType==NonlinearSolverType::BFGS){
        str="  solver type= BFGS";
    }
    else if(_SolverType==NonlinearSolverType::BROYDEN){
        str="  solver type= Broyden";
    }
    else if(_SolverType==NonlinearSolverType::BADBROYDEN){
        str="  solver type= Bad Broyden";
    }
    else if(_SolverType==NonlinearSolverType::NEWTONCG){
        str="  solver type= newton CG";
    }
    else if(_SolverType==NonlinearSolverType::NEWTONGMRES){
        str="  solver type= newton GMRES";
    }
    MessagePrinter::PrintNormalTxt(str);

    snprintf(buff,70,"  max iters=%3d, abs R tol=%13.5e, rel R tol=%13.5e",_MaxIters,_RAbsTol,_RRelTol);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);

    str="  linear solver is: "+_LinearSolverName;
    MessagePrinter::PrintNormalTxt(str);
    
    MessagePrinter::PrintDashLine();
}
