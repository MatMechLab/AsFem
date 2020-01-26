#include "MessagePrint/MessagePrint.h"

void Msg_BCBlockBlockNameInvalid(string blockname)
{
    printf("*** Error: boundname invalid in [%12s]/[bcs]       ***\n",blockname.c_str());
}