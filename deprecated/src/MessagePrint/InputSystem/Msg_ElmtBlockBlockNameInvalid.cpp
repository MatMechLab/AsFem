#include "MessagePrint/MessagePrint.h"


void Msg_ElmtBlockBlockNameInvalid(string blockname)
{
    printf("*** Error: blockname invalid in [%12s]/[elmts]     ***\n",blockname.c_str());
}