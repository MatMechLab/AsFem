#include "FEProblem/FEProblem.h"

void FEProblem::ReadInputFile(){
    inputSystem.ReadInputFile(mesh,dofHandler,
                              bcSystem,icSystem,elmtSystem,materialSystem,
                              qpBlockInfo,
                              linearSolverBlock,
                              nonlinearSolverBlock,
                              solutionSystem,
                              feCtrlInfo);
}