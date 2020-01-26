#include "MessagePrint/MessagePrint.h"

void Msg_Input_JobBlock_InvalidTimeSteppingOption(){
    cout<<"*** Error: invalid time stepping method in input file!     ***"<<endl;
    cout<<"***        timemethod=be[cn,bdf2 ...] is expected !!!      ***"<<endl;
}