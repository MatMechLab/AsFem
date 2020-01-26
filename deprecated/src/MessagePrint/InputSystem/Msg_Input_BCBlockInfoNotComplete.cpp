#include "MessagePrint/MessagePrint.h"

void Msg_Input_BCBlockInfoNotComplete()
{
    cout<<"*** Error: [bcs] block information is not complete!!!      ***"<<endl;
}

//***********************************
void Msg_Input_BCBlockNoDof()
{
    cout<<"***        no 'dof=' found in [bcs] block !!!              ***"<<endl;
}
//***********************************
void Msg_Input_BCBlockNoElmt()
{
    cout<<"***        no 'type=' found in [bcs] block !!!             ***"<<endl;
    
}
//***********************************
void Msg_Input_BCBlockNoBoundary()
{
    cout<<"***        no 'boundary=' found in [bcs] block !!!         ***"<<endl;
}