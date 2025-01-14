//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.08.12
//+++ Purpose: the wrapper for SNES solver (PETSc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "NonlinearSolver/SNESSolver.h"

SNESSolver::SNESSolver(){
    m_Initialized=false;/**< boolean flag for the status of initializing */
    m_MaxIters=50;/**< the maximum iterations */
    m_Iterations=0;/**< the iteration number */
    m_AbsTolR=1.0e-7;/**< the absolute tolerance for residual */
    m_RelTolR=1.0e-9;/**< the relative tolerance for residual */

    m_AbsTolDu=1.0e-4;/**< the absolute tolerance for delta u */
    m_RelTolDu=1.0e-7;/**< the relative tolerance for delta u */

    m_AbsTolE=1.0e-10;/**< the absolute tolerance for energy */
    m_RelTolE=1.0e-15;/**< the relative tolerance for energy */

    m_RNorm0=1.0;/**< the initial norm of residual */
    m_RNorm =1.0;/**< the intermediate or final norm of reisudal */

    m_STol=0.0;

    m_NLSolverName="newton with line search";/**< the nonlinear solver name in SNES */
    m_NLSolverType=NonlinearSolverType::NEWTONLS;
}
void SNESSolver::setFromNonlinearSolverBlock(const NonlinearSolverBlock &nlblock){
    m_MaxIters=nlblock.m_MaxIters;
    m_AbsTolR=nlblock.m_AbsTolR;
    m_RelTolR=nlblock.m_RelTolR;

    m_NLSolverName=nlblock.m_NlSolverTypeName;
    m_NLSolverType=nlblock.m_NlSolverType;

    m_STol=nlblock.m_STol;
}

void SNESSolver::initSolver(LinearSolver &lsolver){
    SNESCreate(PETSC_COMM_WORLD,&m_SNES);

    //**************************************************
    //*** setup KSP
    //**************************************************
    SNESSetKSP(m_SNES,lsolver.getKSPRef());

    //**************************************************
    //*** basic settings for SNES
    //**************************************************
    SNESSetTolerances(m_SNES,m_AbsTolR,m_RelTolR,m_STol,m_MaxIters,-1);
    SNESSetDivergenceTolerance(m_SNES,-1);

    //**************************************************
    //*** for different types of SNES solver
    //**************************************************
    if(m_NLSolverType==NonlinearSolverType::NEWTON||
       m_NLSolverType==NonlinearSolverType::NEWTONLS){
        SNESSetType(m_SNES,SNESNEWTONLS);
        SNESGetLineSearch(m_SNES,&m_SNESLineSearch);
        SNESLineSearchSetType(m_SNESLineSearch,SNESLINESEARCHBT);
        SNESLineSearchSetOrder(m_SNESLineSearch,3);
    }
    else if(m_NLSolverType==NonlinearSolverType::NEWTONAL){
        SNESSetType(m_SNES,SNESNEWTONAL);
    }
    else if(m_NLSolverType==NonlinearSolverType::NEWTONSECANT){
        SNESSetType(m_SNES,SNESNEWTONLS);
        SNESGetLineSearch(m_SNES,&m_SNESLineSearch);
        SNESLineSearchSetType(m_SNESLineSearch,SNESLINESEARCHL2);
    }
    else if(m_NLSolverType==NonlinearSolverType::NEWTONTR){
        SNESSetType(m_SNES,SNESNEWTONTR);
    }
    else if(m_NLSolverType==NonlinearSolverType::BFGS){
        SNESSetType(m_SNES,SNESQN);
    }
    else if(m_NLSolverType==NonlinearSolverType::BROYDEN){
        SNESSetType(m_SNES,SNESQN);
        SNESQNSetType(m_SNES,SNES_QN_BROYDEN);
    }
    else if(m_NLSolverType==NonlinearSolverType::BADBROYDEN){
        SNESSetType(m_SNES,SNESQN);
        SNESQNSetType(m_SNES,SNES_QN_BADBROYDEN);
    }
    else if(m_NLSolverType==NonlinearSolverType::NEWTONCG){
        SNESSetType(m_SNES,SNESNCG);
    }
    else if(m_NLSolverType==NonlinearSolverType::NEWTONGMRES){
        SNESSetType(m_SNES,SNESNGMRES);
    }
    else if(m_NLSolverType==NonlinearSolverType::RICHARDSON){
        SNESSetType(m_SNES,SNESNRICHARDSON);
    }
    else if(m_NLSolverType==NonlinearSolverType::NMS){
        SNESSetType(m_SNES,SNESMS);
        SNESMSSetType(m_SNES,SNESMSEULER);
        PCSetType(lsolver.getPCRef(),PCMG);
    }
    else if(m_NLSolverType==NonlinearSolverType::FAS){
        SNESSetType(m_SNES,SNESFAS);
    }

    SNESSetFromOptions(m_SNES);

    m_Initialized=true;
}

void SNESSolver::releaseMemory(){
    if(m_Initialized){
        SNESDestroy(&m_SNES);
        m_Initialized=false;
    }
}

void SNESSolver::printSolverInfo()const{
    MessagePrinter::printNormalTxt("Nonlinear (SNES) solver information summary:");
    char buff[70];
    string str;

    if(m_NLSolverType==NonlinearSolverType::NEWTON||
       m_NLSolverType==NonlinearSolverType::NEWTONLS){
        str="  Solver type= newton with line search";
    }
    else if(m_NLSolverType==NonlinearSolverType::NEWTONTR){
        str="  Solver type= newton trust region";
    }
    else if(m_NLSolverType==NonlinearSolverType::BFGS){
        str="  Solver type= BFGS";
    }
    else if(m_NLSolverType==NonlinearSolverType::BROYDEN){
        str="  Solver type= Broyden";
    }
    else if(m_NLSolverType==NonlinearSolverType::BADBROYDEN){
        str="  Solver type= Bad Broyden";
    }
    else if(m_NLSolverType==NonlinearSolverType::NEWTONCG){
        str="  Solver type= newton CG";
    }
    else if(m_NLSolverType==NonlinearSolverType::NEWTONGMRES){
        str="  Solver type= newton GMRES";
    }
    MessagePrinter::printNormalTxt(str);

    snprintf(buff,70,"  Max iterations=%3d,",m_MaxIters);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    snprintf(buff,70,"  Abusolute |R| tolerance=%14.5e",m_AbsTolR);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    snprintf(buff,70,"  Relative |R| tolerance=%14.5e",m_RelTolR);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    MessagePrinter::printStars();
    
}