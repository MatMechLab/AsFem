#include "NonLinearSolver/NonLinearSolver.h"

bool NonLinearSolver::Solve(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo){
    switch(GetNonlinearSolverType()){
        case NonLinearSolverType::LINEAR:
            return Linear(mesh,dofHandler,solutionSystem,equationSystem,
                          elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo);
        case NonLinearSolverType::NEWTONRAPHSON:
            return NewtonRaphson(mesh,dofHandler,solutionSystem,equationSystem,
                                 elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo);
        case NonLinearSolverType::MODIFIEDNEWTON:
            return NewtonRaphsonModify(mesh,dofHandler,solutionSystem,equationSystem,
                                       elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo);
        case NonLinearSolverType::NEWTONWITHLINESEARCH:
            return NewtonRaphsonWithLineSearch(mesh,dofHandler,solutionSystem,equationSystem,
                                        elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo);
        case NonLinearSolverType::NEWTONWITHBISECTIONLINESEARCH:
            return NewtonRaphsonWithBisectionLineSearch(mesh,dofHandler,solutionSystem,equationSystem,
                                        elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo);
        case NonLinearSolverType::NEWTONWITHREGULARFALSILINESEARCH:
            return NewtonRaphsonWithRegulaFalsiLineSearch(mesh,dofHandler,solutionSystem,equationSystem,
                                        elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo);
        default:
            assert("*** Error: sorry, unsupported nonlinear solver method!!!");
            break;
    }
    return false;
}