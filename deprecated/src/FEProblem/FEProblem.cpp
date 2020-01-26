#include "FEProblem/FEProblem.h"


FEProblem::FEProblem(int args, char *argv[]){
    JobType=FEJobType::STATIC;
    inputSystem.InitInputSystem(args,argv);
}