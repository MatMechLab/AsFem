#include "FEProblem/FEProblem.h"

void FEProblem::RunFEProblem(){
    StartJob();
    PreRunFEProblem();
    InitFEProblem();
    RunFEAnalysis();
    Finalize();
}