#include "FEProblem/FEProblem.h"

void FEProblem::RunFEAnalysis(){
    PrintJobInfo();
    if(JobType==FEJobType::STATIC){
        RunStaticAnalysis();
    }
    else{
        RunTransientAnalysis();
    }
}