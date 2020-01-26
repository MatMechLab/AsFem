#include "NonLinearSolver/NonLinearSolver.h"

void NonLinearSolver::InitNonlinearSolver(NonLinearSolverBlockInfo &nonlinearSolverBlockInfo,
                             LinearSolverBlockInfo &linearSolverBlockInfo){
   
    _SolverType=nonlinearSolverBlockInfo._SolverType;

    _IsDebug=false;
    _ConvergenceCriterion=nonlinearSolverBlockInfo._Criterion;
                            //0-->use the residual or energy
                            //1-->only use the residual
                            //2-->only use the energy
                            //3-->use the s tol of line search algorithm

    _CurrentIters=0;
    _MaxIters=nonlinearSolverBlockInfo._MaxIters;
    _MaxLineSearch=nonlinearSolverBlockInfo._MaxLineSearch;

    _RAbsTol=nonlinearSolverBlockInfo._Rabstol;
    _RRelTol=nonlinearSolverBlockInfo._Rreltol;

    _EAbsTol=nonlinearSolverBlockInfo._Eabstol;
    _ERelTol=nonlinearSolverBlockInfo._Ereltol;

    _Stol=nonlinearSolverBlockInfo._Stol;

    _Rnorm0=1.0;_Rnorm=1.0;
    _Enorm0=1.0;_Enorm=1.0;

    // for linear solver
    _linearSolver.InitLinearSolver(linearSolverBlockInfo);
}