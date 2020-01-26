#include "FEProblem/FEProblem.h"

void FEProblem::RunTransientAnalysis(){
    timeStepping.RunTimeStepping(mesh,dofHandler,
        solutionSystem,equationSystem,
        elmtSystem,materialSystem,
        bcSystem,icSystem,
        feSystem,nonlinearSolver,feCtrlInfo,
        outputSystem);
}