#include "FEProblem/FEProblem.h"

void FEProblem::PreRunFEProblem(){
    
    cout<<"**************************************************************"<<endl;
    cout<<"*** Start to read input file ...                           ***"<<endl;
    _TimerStartOfInput=chrono::high_resolution_clock::now();
    ReadInputFile();
    _TimerEndOfInput=chrono::high_resolution_clock::now();
    _DurationOfInput=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfInput-_TimerStartOfInput).count()/1.0e6;
    printf("*** Read input file finished!          ===>[%12.5f s]***\n",_DurationOfInput);
    cout<<"***--------------------------------------------------------***"<<endl;
    
    //**************************************************************
    //*** Check the input file settings is correct or not
    //**************************************************************
    cout<<"*** Start to check input file settings ...                 ***"<<endl;
    _TimerStartOfInputCheck=chrono::high_resolution_clock::now();
    CheckInputFileInfo();
    _TimerEndOfInputChekc=chrono::high_resolution_clock::now();
    _DurationOfInputCheck=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfInputChekc-_TimerStartOfInputCheck).count()/1.0e6;
    printf("*** Check input file finished!         ===>[%12.5f s]***\n",_DurationOfInputCheck);
    cout<<"***--------------------------------------------------------***"<<endl;

}