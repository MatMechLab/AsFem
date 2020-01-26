#include "NonLinearSolver/NonLinearSolver.h"

bool NonLinearSolver::ConvergenceCheck(){
    //0-->use the residual or energy
    //1-->only use the residual
    //2-->only use the energy
    //3-->use both the residual and energy
    switch (GetConvergenceCriterion())
    {
    case 0:
        if((_Rnorm<_RAbsTol||_Rnorm<_Rnorm0*_RRelTol)||
            (_Enorm<_EAbsTol||_Enorm<_Enorm0*_ERelTol)){
            return true;
        }
        else{
            return false;
        }
        break;
    case 1:
        if(_Rnorm<_RAbsTol||_Rnorm<_Rnorm0*_RRelTol){
            return true;
        }
        else{
            return false;
        }
        break;
    case 2:
        if(_Enorm<_EAbsTol||_Enorm<_Enorm0*_ERelTol){
            return true;
        }
        else{
            return false;
        }
        break;
    case 3:
        if((_Rnorm<_RAbsTol||_Rnorm<_Rnorm0*_RRelTol)&&
            (_Enorm<_EAbsTol||_Enorm<_Enorm0*_ERelTol)){
            return true;
        }
        else{
            return false;
        }
        break;
    default:
        return false;
    }
    return false;
}