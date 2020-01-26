#include "LinearSolver/LinearSolver.h"

void LinearSolver::PrintSolverInfo() const{
    switch (GetSolverType())
    {
        case LinearSolverType::SPARSELU:
            printf("***   linear solver='SparseLU'                             ***\n");
            break;
        case LinearSolverType::SPARSEQR:
            printf("***   linear solver='SparseQR'                             ***\n");
            break;
        case LinearSolverType::PARDISOLU:
            printf("***   linear solver='PARDISO'                              ***\n");
            break;
        case LinearSolverType::BICG:
            printf("***   linear solver='BICG'                                 ***\n");
            printf("***   max iters=%8d, tolerance=%14.6e         ***\n",_MaxIters,_Tolerance);
            break;
        case LinearSolverType::CG:
            printf("***   linear solver='CG'                                   ***\n");
            printf("***   max iters=%8d, tolerance=%14.6e         ***\n",_MaxIters,_Tolerance);
            break;
        case LinearSolverType::GMRES:
            printf("***   linear solver='GMRES', restart=%8d              ***\n",_Restarts);
            printf("***   max iters=%8d, tolerance=%14.6e         ***\n",_MaxIters,_Tolerance);
            break;
        case LinearSolverType::UMFPACK:
            printf("***   linear solver='UMFPACK of SuiteSparse'               ***\n");
            printf("***   max iters=%8d, tolerance=%14.6e         ***\n",_MaxIters,_Tolerance);
            break;
        default:
            break;
    }
}