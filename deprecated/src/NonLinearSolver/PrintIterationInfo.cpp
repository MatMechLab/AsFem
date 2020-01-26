#include "NonLinearSolver/NonLinearSolver.h"

void NonLinearSolver::PrintDepIterationInfo() const{
    if(GetNonlinearSolverType()==NonLinearSolverType::NEWTONWITHLINESEARCH){
        printf("***    Iters=%3d, |R|=%12.6e, |E|=%12.6e       ***\n",_CurrentIters,_Rnorm,_Enorm);
        printf("***                 S=%12.5e, Eta=%12.5e       ***\n",_s,_eta);
    }
    else{
        printf("***    Iters=%3d, |R|=%12.6e, |E|=%12.6e       ***\n",_CurrentIters,_Rnorm,_Enorm);
    }
}
//*********************************
void NonLinearSolver::PrintIterationInfo() const{
    if(GetNonlinearSolverType()==NonLinearSolverType::NEWTONWITHLINESEARCH){
        printf("***    Iters=%3d:|R0|=%12.6e, |R|=%12.6e       ***\n",_CurrentIters,_Rnorm0,_Rnorm);
        printf("***              |E0|=%12.6e, |E|=%12.6e       ***\n",_Enorm0,_Enorm);
        printf("***                 S=%12.6e, Eta=%12.6e       ***\n",_s,_eta);
    }
    else{
        printf("***    Iters=%3d:|R0|=%12.6e, |R|=%12.6e       ***\n",_CurrentIters,_Rnorm0,_Rnorm);
        printf("***              |E0|=%12.6e, |E|=%12.6e       ***\n",_Enorm0,_Enorm);
    }
}