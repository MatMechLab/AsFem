#include "MessagePrint/MessagePrint.h"

void Msg_Input_ICBlockNoValueWarning(string blockname)
{
    printf("*** Warning: no values in [%18s]/[ics]!!     ***\n",blockname.c_str());
    cout<<"***           default initial value will be used !!!       ***"<<endl;
}