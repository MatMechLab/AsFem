#include "MessagePrint/MessagePrint.h"

void Msg_ICBlockBlockNameInvalid(string blockname)
{
    printf("*** Error: blockname invalid in [%12s]/[ics]       ***\n",blockname.c_str());
}