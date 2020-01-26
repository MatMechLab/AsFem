#include "MessagePrint/MessagePrint.h"

void Msg_Input_ElmtBlockInfoNotComplete()
{
    cout<<"*** Error: [elmts] block information is not complete!!     ***"<<endl;
    
}
//***********************************
void Msg_Input_ElmtBlockNoType()
{
    cout<<"***        no 'type=' found in [elmts] block !!!           ***"<<endl;
}
//***********************************
void Msg_Input_ElmtBlockNoDofs()
{
    cout<<"***        no 'dofs=' found in [elmts] block !!!           ***"<<endl;
}