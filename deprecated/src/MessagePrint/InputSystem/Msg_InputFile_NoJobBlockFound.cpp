#include "MessagePrint/MessagePrint.h"

void Msg_InputFile_NoJobBlockFound()
{
    cout<<"*** Error: no [job] block found in your input file!!!      ***"<<endl;
    cout<<"***        [job]/[end] block must be given for AsFem!!     ***"<<endl;
}