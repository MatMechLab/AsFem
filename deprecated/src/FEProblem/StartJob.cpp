#include "FEProblem/FEProblem.h"

void FEProblem::StartJob(){
    _JobStartTime=chrono::system_clock::now();
    auto in_time_t=chrono::system_clock::to_time_t(_JobStartTime);
    stringstream ss;

    ss<<put_time(localtime(&in_time_t),"%Y-%m-%d %X");
    _JobStartTimeContext=ss.str();
    
    printf("*** AsFem job starts at: %25s         ***\n",_JobStartTimeContext.c_str());
}