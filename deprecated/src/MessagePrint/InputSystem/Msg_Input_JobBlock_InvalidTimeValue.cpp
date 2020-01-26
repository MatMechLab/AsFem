#include "MessagePrint/MessagePrint.h"

void Msg_Input_JobBlock_InvalidTimeValue(){
    cout<<"*** Error: invalid endtime value in input file!!!          ***"<<endl;
    cout<<"***        endtime=10^-10~10^20 is expected in [job]!!     ***"<<endl;
}