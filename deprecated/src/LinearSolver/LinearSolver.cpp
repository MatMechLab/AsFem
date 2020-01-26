#include "LinearSolver/LinearSolver.h"


LinearSolver::LinearSolver(){
    _SolverType=LinearSolverType::SPARSELU;
    _MaxIters=20000;
    _Restarts=500;
    _Tolerance=1.0e-10;
    _IsInit=false;
}