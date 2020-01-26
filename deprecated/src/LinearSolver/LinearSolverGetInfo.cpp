#include "LinearSolver/LinearSolver.h"

int LinearSolver::GetIters() const
{
    switch (GetSolverType())
    {
        case LinearSolverType::BICG:
            return _BICGSolver.iterations();
        case LinearSolverType::CG:
            return _CGSolver.iterations();
        default:
            return 1;
    }
}
//**********************************
double LinearSolver::GetTolerance() const{
    switch (GetSolverType())
    {
        case LinearSolverType::BICG:
            return _BICGSolver.tolerance();
        case LinearSolverType::CG:
            return _CGSolver.tolerance();
        default:
            return 1.0e-16;
    }
}