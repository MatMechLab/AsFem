#include "MessagePrint/MessagePrint.h"

void Msg_Input_MateBlockInfoNotComplete()
{
    cout<<"*** Error: [mate] block information is not complete!!!     ***"<<endl;
}

//***********************************
void Msg_Input_MateBlockNoElmt()
{
    cout<<"***        no 'type=' found in [mate] block !!!            ***"<<endl;

}
//***********************************
void Msg_Input_MateBlockNoParams()
{
    cout<<"***        no 'params=' found in [mate] block !!!          ***"<<endl;
}