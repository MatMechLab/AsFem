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

#include "NonlinearSolver/NonlinearSolver.h"

NonlinearSolver::NonlinearSolver(){
    _Rnorm0=1.0;_Rnorm=1.0;
    _Enorm0=1.0;_Enorm=1.0;
    _RAbsTol=1.0e-8;_RRelTol=1.0e-10;
    _EAbsTol=1.0e-19;_ERelTol=1.0e-20;
    _MaxIters=20;_Iters=0;
    _STol=1.0e-16;
    _SolverType=NonlinearSolverType::NEWTONLS;
    _SolverName="newton with line search";
    _PCTypeName="lu";
}

void NonlinearSolver::SetOptionsFromNonlinearSolverBlock(NonlinearSolverBlock &nonlinearsolverblock){
    _SolverType=nonlinearsolverblock._SolverType;
    _MaxIters=nonlinearsolverblock._MaxIters;
    _RAbsTol=nonlinearsolverblock._RAbsTol;
    _RRelTol=nonlinearsolverblock._RRelTol;
    _STol=nonlinearsolverblock._STol;

    _SolverType=nonlinearsolverblock._SolverType;
    _SolverName=nonlinearsolverblock._SolverTypeName;
    _PCTypeName=nonlinearsolverblock._PCTypeName;
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
    KSPGMRESSetRestart(_ksp,1300);
    KSPGetPC(_ksp,&_pc);
    #ifdef HASMUMPS
    PCSetType(_pc,PCLU);
    KSPSetType(_ksp,KSPPREONLY);
    PCFactorSetMatSolverType(_pc,MATSOLVERMUMPS);
    #endif

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

//***************************************************
void NonlinearSolver::ReleaseMem(){
    SNESDestroy(&_snes);
}

//****************************************************
void NonlinearSolver::PrintInfo()const{
    MessagePrinter::PrintNormalTxt(" Nonlinear solver information summary:");
    char buff[70];
    string str;

    if(_SolverType==NonlinearSolverType::NEWTON||_SolverType==NonlinearSolverType::NEWTONLS){
        str="   solver type= newton with line search";
    }
    else if(_SolverType==NonlinearSolverType::NEWTONTR){
        str="   solver type= newton trust region";
    }
    else if(_SolverType==NonlinearSolverType::BFGS){
        str="   solver type= BFGS";
    }
    else if(_SolverType==NonlinearSolverType::BROYDEN){
        str="   solver type= Broyden";
    }
    else if(_SolverType==NonlinearSolverType::BADBROYDEN){
        str="   solver type= Bad Broyden";
    }
    else if(_SolverType==NonlinearSolverType::NEWTONCG){
        str="   solver type= newton CG";
    }
    else if(_SolverType==NonlinearSolverType::NEWTONGMRES){
        str="   solver type= newton GMRES";
    }
    MessagePrinter::PrintNormalTxt(str);

    snprintf(buff,70,"   max iters=%3d, abs R tol=%13.5e, rel R tol=%13.5e",_MaxIters,_RAbsTol,_RRelTol);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    
    MessagePrinter::PrintDashLine();
}