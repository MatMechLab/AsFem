#include "LinearSolver/LinearSolver.h"


void LinearSolver::InitLinearSolver(LinearSolverBlockInfo &linearSolverBlockInfo){
    _SolverType=linearSolverBlockInfo._SolverType;
    _MaxIters=linearSolverBlockInfo._MaxIters;
    _Tolerance=linearSolverBlockInfo._Tol;
    _Restarts=linearSolverBlockInfo._Restart;
    _IsInit=false;
}