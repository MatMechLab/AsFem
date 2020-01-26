#include "MessagePrint/MessagePrint.h"

void Msg_Input_ICBlockInfoNotComplete()
{
    cout<<"*** Error: [ics] block information is not complete!!!      ***"<<endl;
}

//***********************************
void Msg_Input_ICBlockNoDof()
{
    cout<<"***        no 'dof=' found in [ics] block !!!              ***"<<endl;
}
//***********************************
void Msg_Input_ICBlockNoElmt()
{
    cout<<"***        no 'type=' found in [ics] block !!!             ***"<<endl;
    
}
//***********************************
void Msg_Input_ICBlockNoBlock()
{
    cout<<"***        no 'block=' found in [ics] block !!!            ***"<<endl;
}