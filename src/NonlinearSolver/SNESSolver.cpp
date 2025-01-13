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
    m_initialized=false;/**< boolean flag for the status of initializing */
    m_maxiters=50;/**< the maximum iterations */
    m_iterations=0;/**< the iteration number */
    m_abstol_r=1.0e-7;/**< the absolute tolerance for residual */
    m_reltol_r=1.0e-9;/**< the relative tolerance for residual */

    m_abstol_du=1.0e-4;/**< the absolute tolerance for delta u */
    m_reltol_du=1.0e-7;/**< the relative tolerance for delta u */

    m_abstol_e=1.0e-10;/**< the absolute tolerance for energy */
    m_reltol_e=1.0e-15;/**< the relative tolerance for energy */

    m_rnorm0=1.0;/**< the initial norm of residual */
    m_rnorm =1.0;/**< the intermediate or final norm of reisudal */

    m_s_tol=0.0;

    m_nlsolvername="newton with line search";/**< the nonlinear solver name in SNES */
    m_nlsolvertype=NonlinearSolverType::NEWTONLS;
}
void SNESSolver::setFromNonlinearSolverBlock(const NonlinearSolverBlock &nlblock){
    m_maxiters=nlblock.m_MaxIters;
    m_abstol_r=nlblock.m_AbsTolR;
    m_reltol_r=nlblock.m_RelTolR;

    m_nlsolvername=nlblock.m_NlSolverTypeName;
    m_nlsolvertype=nlblock.m_NlSolverType;

    m_s_tol=nlblock.m_STol;
}

void SNESSolver::initSolver(LinearSolver &lsolver){
    SNESCreate(PETSC_COMM_WORLD,&m_snes);

    //**************************************************
    //*** setup KSP
    //**************************************************
    SNESSetKSP(m_snes,lsolver.getKSPRef());

    //**************************************************
    //*** basic settings for SNES
    //**************************************************
    SNESSetTolerances(m_snes,m_abstol_r,m_reltol_r,m_s_tol,m_maxiters,-1);
    SNESSetDivergenceTolerance(m_snes,-1);

    //**************************************************
    //*** for different types of SNES solver
    //**************************************************
    if(m_nlsolvertype==NonlinearSolverType::NEWTON||
       m_nlsolvertype==NonlinearSolverType::NEWTONLS){
        SNESSetType(m_snes,SNESNEWTONLS);
        SNESGetLineSearch(m_snes,&m_sneslinesearch);
        SNESLineSearchSetType(m_sneslinesearch,SNESLINESEARCHBT);
        SNESLineSearchSetOrder(m_sneslinesearch,3);
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONAL){
        SNESSetType(m_snes,SNESNEWTONAL);
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONSECANT){
        SNESSetType(m_snes,SNESNEWTONLS);
        SNESGetLineSearch(m_snes,&m_sneslinesearch);
        SNESLineSearchSetType(m_sneslinesearch,SNESLINESEARCHL2);
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONTR){
        SNESSetType(m_snes,SNESNEWTONTR);
    }
    else if(m_nlsolvertype==NonlinearSolverType::BFGS){
        SNESSetType(m_snes,SNESQN);
    }
    else if(m_nlsolvertype==NonlinearSolverType::BROYDEN){
        SNESSetType(m_snes,SNESQN);
        SNESQNSetType(m_snes,SNES_QN_BROYDEN);
    }
    else if(m_nlsolvertype==NonlinearSolverType::BADBROYDEN){
        SNESSetType(m_snes,SNESQN);
        SNESQNSetType(m_snes,SNES_QN_BADBROYDEN);
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONCG){
        SNESSetType(m_snes,SNESNCG);
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONGMRES){
        SNESSetType(m_snes,SNESNGMRES);
    }
    else if(m_nlsolvertype==NonlinearSolverType::RICHARDSON){
        SNESSetType(m_snes,SNESNRICHARDSON);
    }
    else if(m_nlsolvertype==NonlinearSolverType::NMS){
        SNESSetType(m_snes,SNESMS);
        SNESMSSetType(m_snes,SNESMSEULER);
        PCSetType(lsolver.getPCRef(),PCMG);
    }
    else if(m_nlsolvertype==NonlinearSolverType::FAS){
        SNESSetType(m_snes,SNESFAS);
    }

    SNESSetFromOptions(m_snes);

    m_initialized=true;
}

void SNESSolver::releaseMemory(){
    if(m_initialized){
        SNESDestroy(&m_snes);
        m_initialized=false;
    }
}

void SNESSolver::printSolverInfo()const{
    MessagePrinter::printNormalTxt("Nonlinear (SNES) solver information summary:");
    char buff[70];
    string str;

    if(m_nlsolvertype==NonlinearSolverType::NEWTON||
       m_nlsolvertype==NonlinearSolverType::NEWTONLS){
        str="  Solver type= newton with line search";
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONTR){
        str="  Solver type= newton trust region";
    }
    else if(m_nlsolvertype==NonlinearSolverType::BFGS){
        str="  Solver type= BFGS";
    }
    else if(m_nlsolvertype==NonlinearSolverType::BROYDEN){
        str="  Solver type= Broyden";
    }
    else if(m_nlsolvertype==NonlinearSolverType::BADBROYDEN){
        str="  Solver type= Bad Broyden";
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONCG){
        str="  Solver type= newton CG";
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONGMRES){
        str="  Solver type= newton GMRES";
    }
    MessagePrinter::printNormalTxt(str);

    snprintf(buff,70,"  Max iterations=%3d,",m_maxiters);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    snprintf(buff,70,"  Abusolute |R| tolerance=%14.5e",m_abstol_r);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    snprintf(buff,70,"  Relative |R| tolerance=%14.5e",m_reltol_r);
    str=buff;
    MessagePrinter::printNormalTxt(str);
    MessagePrinter::printStars();
    
}