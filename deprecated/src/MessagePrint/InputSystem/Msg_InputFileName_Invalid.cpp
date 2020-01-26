#include "MessagePrint/MessagePrint.h"

void Msg_InputFileName_Invalid(string filename)
{
    printf("*** Error: can't open  %-30s      ***\n",filename.c_str());
}