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
//+++ Date   : 2025.01.11
//+++ Purpose: This class offers the API for the KSP linear solver
//+++          from PETSc, one should have PETSc installed
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "LinearSolver/KSPSolver.h"

KSPSolver::KSPSolver() {
    m_Iterations=0;
    m_MaxIterations=10000;
    m_GMRESRestartNumber=1500;
    m_Tolerance=1.0e-25;

    m_KSPSolverTypeName="gmres";
    string m_PCTypeName="bjacobi";
}
void KSPSolver::setDefaultParams(){
    m_MaxIterations=10000;
    m_GMRESRestartNumber=1500;
    m_Tolerance=1.0e-16;

    m_KSPSolverTypeName="gmres";
    string m_PCTypeName="bjacobi";
}
void KSPSolver::init() {
    KSPCreate(PETSC_COMM_WORLD, &m_KSP);
    KSPGetPC(m_KSP, &m_PC);

    KSPGMRESSetRestart(m_KSP,m_GMRESRestartNumber);
    KSPSetTolerances(m_KSP,PETSC_DEFAULT,m_Tolerance,PETSC_DEFAULT,m_MaxIterations);

    // setup the preconditioner
    if (m_PCTypeName == "jacobi") {
        PCSetType(m_PC,PCJACOBI);
    }
    else if (m_PCTypeName == "bjacobi") {
        PCSetType(m_PC,PCBJACOBI);
    }
    else if (m_PCTypeName == "sor") {
        PCSetType(m_PC,PCSOR);
    }
    else if (m_PCTypeName == "eisenstat") {
        PCSetType(m_PC,PCEISENSTAT);
    }
    else if (m_PCTypeName == "icc") {
        PCSetType(m_PC,PCICC);
    }
    else if (m_PCTypeName == "ilu") {
        PCSetType(m_PC,PCILU);
    }
    else if (m_PCTypeName == "asm") {
        PCSetType(m_PC,PCASM);
    }
    else if (m_PCTypeName == "gasm") {
        PCSetType(m_PC,PCGASM);
    }
    else if (m_PCTypeName == "gamg") {
        PCSetType(m_PC,PCGAMG);
    }
    else if (m_PCTypeName == "bddc") {
        PCSetType(m_PC,PCBDDC);
    }
    else if (m_PCTypeName == "ksp") {
        PCSetType(m_PC,PCKSP);
    }
    else if (m_PCTypeName == "composite") {
        PCSetType(m_PC,PCCOMPOSITE);
    }
    else if (m_PCTypeName == "lu") {
        PCSetType(m_PC,PCLU);
        PCFactorSetZeroPivot(m_PC,1.0e-15);
    }
    else if (m_PCTypeName == "cholesky") {
        PCSetType(m_PC,PCCHOLESKY);
    }
    else if (m_PCTypeName == "none") {
        PCSetType(m_PC,PCNONE);
    }
    else if (m_PCTypeName == "shell") {
        PCSetType(m_PC,PCSHELL);
    }
    else {
        MessagePrinter::printTxt("preconditoner="+m_PCTypeName+" is invalid, please check your input file");
        MessagePrinter::exitAsFem();
    }

    if (m_KSPSolverTypeName == "gmres") {
        KSPSetType(m_KSP,KSPGMRES);
    }
    else if (m_KSPSolverTypeName == "fgmres") {
        KSPSetType(m_KSP,KSPFGMRES);
    }
    else if (m_KSPSolverTypeName == "lgmres") {
        KSPSetType(m_KSP,KSPLGMRES);
    }
    else if (m_KSPSolverTypeName == "dgmres") {
        KSPSetType(m_KSP,KSPDGMRES);
    }
    else if (m_KSPSolverTypeName == "pgmres") {
        KSPSetType(m_KSP,KSPPGMRES);
    }
    else if (m_KSPSolverTypeName == "cg") {
        KSPSetType(m_KSP,KSPCG);
    }
    else if (m_KSPSolverTypeName == "chebyshev") {
        KSPSetType(m_KSP,KSPCHEBYSHEV);
    }
    else if (m_KSPSolverTypeName == "richardson") {
        KSPSetType(m_KSP,KSPRICHARDSON);
    }
    else if (m_KSPSolverTypeName == "bcgs") {
        KSPSetType(m_KSP,KSPBCGS);//the BiCGStab (Stabilized version of Biconjugate Gradient) method
    }
    else if (m_KSPSolverTypeName == "bicg") {
        KSPSetType(m_KSP,KSPBICG);
    }
    else if (m_KSPSolverTypeName == "mumps") {
        KSPSetType(m_KSP,KSPNONE);
        PCSetType(m_PC,PCLU);
        PCFactorSetMatSolverType(m_PC,MATSOLVERMUMPS);
    }
    else if (m_KSPSolverTypeName == "superlu") {
        KSPSetType(m_KSP,KSPNONE);
        PCSetType(m_PC,PCLU);
        PCFactorSetMatSolverType(m_PC,MATSOLVERSUPERLU_DIST);
    }
    KSPSetFromOptions(m_KSP);
    PCSetFromOptions(m_PC);
    KSPSetUp(m_KSP);
}
bool KSPSolver::solve(SparseMatrix &A,Vector &b,Vector &x) {
    m_Iterations=0;
    KSPSetOperators(m_KSP,A.getCopy(),A.getCopy());
    KSPSolve(m_KSP,b.getVectorCopy(),x.getVectorRef());
    KSPGetIterationNumber(m_KSP,&m_Iterations);
    return true;
}
void KSPSolver::printKSPSolverInfo()const {
    MessagePrinter::printTxt("KSP solver info:");
    MessagePrinter::printTxt("  Max iterations=\'"+to_string(m_MaxIterations)+"\'");
    MessagePrinter::printTxt("  Restarts=\'"+to_string(m_GMRESRestartNumber)+"\'");
    MessagePrinter::printTxt("  Solver=\'"+m_KSPSolverTypeName+"\'");
    MessagePrinter::printTxt("  PC type=\'"+m_PCTypeName+"\'");
    char buff[65];
    snprintf(buff,65,"  Tolerance=%14.6e",m_Tolerance);
    MessagePrinter::printTxt(buff);
    MessagePrinter::printStars();
}