#include "MessagePrint/MessagePrint.h"

void Msg_Input_MateBlockNoParamsWarning(string blockname)
{
    printf("*** Warning: no params in [%18s]/[mates]     ***\n",blockname.c_str());
    cout<<"***           default material value will be used !!!      ***"<<endl;
}