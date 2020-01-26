#include "FEProblem/FEProblem.h"

void FEProblem::Finalize(){
    _JobEndTime=chrono::system_clock::now();
    auto in_time_t=chrono::system_clock::to_time_t(_JobEndTime);
    stringstream ss;
    ss<<put_time(localtime(&in_time_t),"%Y-%m-%d %X");
    _JobEndTimeContext=ss.str();
    printf("**************************************************************\n");
    printf("*** AsFem job finish at: %25s         ***\n",_JobEndTimeContext.c_str());
    printf("**************************************************************\n");
    printf("*** Thank you for the use !                                ***\n");
    printf("**************************************************************\n");
}
