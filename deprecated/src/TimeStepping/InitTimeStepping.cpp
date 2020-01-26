#include "TimeStepping/TimeStepping.h"

void TimeStepping::InitTimeStepping(){
    ctan[0]=1.0;ctan[1]=1.0;
    _IsAdaptive=false;
    _Method=TimeSteppingType::BACKWARDEULER;
    _TotalTime=0.0;_dt0=1.0e-6;_dt=1.0e-6;
    _GrowthFactor=1.2;_CutbackFactor=0.75;
    _TotalStep=0;_CurrentStep=0;
}