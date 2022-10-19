//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
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

    m_linearsolvername="gmres";/**< the string name of the linear solver in SNES*/
    m_nlsolvername="newton with line search";/**< the nonlinear solver name in SNES */
    m_pcname="lu";/**< the preconditioner name of current SNES solver */
    m_nlsolvertype=NonlinearSolverType::NEWTONLS;
}
void SNESSolver::setFromNonlinearSolverBlock(const NonlinearSolverBlock &nlblock){
    m_maxiters=nlblock.m_maxiters;
    m_abstol_r=nlblock.m_abstol_r;
    m_reltol_r=nlblock.m_reltol_r;

    m_nlsolvername=nlblock.m_nlsolvertypename;
    m_nlsolvertype=nlblock.m_nlsolvertype;

    m_linearsolvername=nlblock.m_linearsolvername;

    m_s_tol=nlblock.m_s_tol;

    m_pcname=nlblock.m_pctypename;
}

void SNESSolver::initSolver(){
    SNESCreate(PETSC_COMM_WORLD,&m_snes);

    //**************************************************
    //*** init KSP
    //**************************************************
    SNESGetKSP(m_snes,&m_ksp);
    KSPGMRESSetRestart(m_ksp,2500);
    KSPGetPC(m_ksp,&m_pc);


    //**************************************************
    //*** setup the preconditioner
    //**************************************************
    if(m_pcname=="lu"){
        PCSetType(m_pc,PCLU);
    }
    else if(m_pcname=="ilu"){
        PCSetType(m_pc,PCILU);
    }
    else if(m_pcname=="jacobi"){
        PCSetType(m_pc,PCJACOBI);
    }
    else if(m_pcname=="bjacobi"){
        PCSetType(m_pc,PCBJACOBI);
    }
    else if(m_pcname=="sor"){
        PCSetType(m_pc,PCSOR);
    }
    else if(m_pcname=="icc"){
        PCSetType(m_pc,PCICC);
    }
    else if(m_pcname=="asm"){
        PCSetType(m_pc,PCASM);
    }
    else if(m_pcname=="gasm"){
        PCSetType(m_pc,PCGASM);
    }
    else if(m_pcname=="gamg"){
        PCSetType(m_pc,PCGAMG);
    }
    else if(m_pcname=="ksp"){
        PCSetType(m_pc,PCKSP);
    }
    else if(m_pcname=="cholesky"){
        PCSetType(m_pc,PCCHOLESKY);
    }
    else if(m_pcname=="none"){
        PCSetType(m_pc,PCNONE);// no preconditioner
    }
    else{
        MessagePrinter::printErrorTxt("unsupported preconditioner("+m_pcname+") in SNESSolver");
        MessagePrinter::exitAsFem();
    }

    //**************************************************
    //*** setup the linear solver
    //**************************************************
    if(m_linearsolvername=="default"){
        PCSetType(m_pc,PCLU);
    }
    else if(m_linearsolvername=="gmres"){
        KSPSetType(m_ksp,KSPGMRES);
    }
    else if(m_linearsolvername=="fgmres"){
        KSPSetType(m_ksp,KSPFGMRES);
    }
    else if(m_linearsolvername=="cg"){
        KSPSetType(m_ksp,KSPCG);
    }
    else if(m_linearsolvername=="bicg"){
        KSPSetType(m_ksp,KSPBICG);
    }
    else if(m_linearsolvername=="richardson"){
        KSPSetType(m_ksp,KSPRICHARDSON);
    }
    else if(m_linearsolvername=="mumps"){
        KSPSetType(m_ksp,KSPPREONLY);
        PCSetType(m_pc,PCLU);
        PCFactorSetMatSolverType(m_pc,MATSOLVERMUMPS);
    }
    else if(m_linearsolvername=="superlu"){
        KSPSetType(m_ksp,KSPPREONLY);
        PCSetType(m_pc,PCLU);
        PCFactorSetMatSolverType(m_pc,MATSOLVERSUPERLU_DIST);
    }


    PCFactorSetReuseOrdering(m_pc,PETSC_TRUE);

    //*** allow user setting ksp from command line
    //**************************************************
    KSPSetFromOptions(m_ksp);
    PCSetFromOptions(m_pc);

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
        PCSetType(m_pc,PCMG);
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
        str="  solver type= newton with line search";
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONTR){
        str="  solver type= newton trust region";
    }
    else if(m_nlsolvertype==NonlinearSolverType::BFGS){
        str="  solver type= BFGS";
    }
    else if(m_nlsolvertype==NonlinearSolverType::BROYDEN){
        str="  solver type= Broyden";
    }
    else if(m_nlsolvertype==NonlinearSolverType::BADBROYDEN){
        str="  solver type= Bad Broyden";
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONCG){
        str="  solver type= newton CG";
    }
    else if(m_nlsolvertype==NonlinearSolverType::NEWTONGMRES){
        str="  solver type= newton GMRES";
    }
    MessagePrinter::printNormalTxt(str);

    snprintf(buff,70,"  max iters=%3d, abs R tol=%13.5e, rel R tol=%13.5e",m_maxiters,m_abstol_r,m_reltol_r);
    str=buff;
    MessagePrinter::printNormalTxt(str);

    str="  linear solver is: "+m_linearsolvername;
    str+=", preconditioner= "+m_pcname;
    MessagePrinter::printNormalTxt(str);
    MessagePrinter::printStars();
    
}