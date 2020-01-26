#include "MessagePrint/MessagePrint.h"

void Msg_Input_JobBlock_InvalidDtValue(){
    cout<<"*** Error: invalid dt value in input file!!!               ***"<<endl;
    cout<<"***        dt=10^-10~10^20 is expected in [job]!!          ***"<<endl;
}