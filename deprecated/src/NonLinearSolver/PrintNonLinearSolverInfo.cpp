#include "NonLinearSolver/NonLinearSolver.h"

void NonLinearSolver::PrintNonLinearSolverInfo() const{
    
    printf("*** Nonlinear solver information summary:                  ***\n");
    if(GetNonlinearSolverType()==NonLinearSolverType::LINEAR){
        printf("***   solver type='linear'                                 ***\n");
    }
    else if(GetNonlinearSolverType()==NonLinearSolverType::NEWTONRAPHSON){
        printf("***   solver type='newton-raphson'                         ***\n");
        printf("***   maximum iterations = %3d                             ***\n",GetMaxIters());
    }
    else if(GetNonlinearSolverType()==NonLinearSolverType::NEWTONWITHLINESEARCH){
        printf("***   solver type='newton-raphson with line search'        ***\n");
        printf("***   maximum iterations = %3d, s error=%14.6e     ***\n",GetMaxIters(),_Stol);
    }
    else if(GetNonlinearSolverType()==NonLinearSolverType::NEWTONWITHBISECTIONLINESEARCH){
        printf("***   solver type='newton-raphson with bisect line search' ***\n");
        printf("***   maximum iterations = %3d, s error=%14.6e     ***\n",GetMaxIters(),_Stol);
    }
    else if(GetNonlinearSolverType()==NonLinearSolverType::NEWTONWITHREGULARFALSILINESEARCH){
        printf("***   solver type='newton-raphson with regula line search' ***\n");
        printf("***   maximum iterations = %3d, s error=%14.6e     ***\n",GetMaxIters(),_Stol);
    }
    else if(GetNonlinearSolverType()==NonLinearSolverType::MODIFIEDNEWTON){
        printf("***   solver type='modified newton-raphson'                ***\n");
        printf("***   maximum iterations = %3d                             ***\n",GetMaxIters());
    }
    
    printf("***   |R| rel error=%10.5e,abs error=%10.5e      ***\n",GetResidualRelativeError(),
                                                                   GetResidualAbsoluteError());
    printf("***   |E| rel error=%10.5e,abs error=%10.5e      ***\n",GetEnergyRelativeError(),
                                                                   GetEnergyAbsoluteError());
                                                                   
    if(GetConvergenceCriterion()==0){
        printf("***   convergence criterion: use residual or energy        ***\n");
    }
    else if(GetConvergenceCriterion()==1){
        printf("***   convergence criterion: only use residual             ***\n");
    }
    else if(GetConvergenceCriterion()==2){
        printf("***   convergence criterion: only use energy               ***\n");
    }
    else if(GetConvergenceCriterion()==3){
        printf("***   convergence criterion: use residual and energy       ***\n");
    }
    _linearSolver.PrintSolverInfo();
    printf("***--------------------------------------------------------***\n");
}