#include "MessagePrint/MessagePrint.h"


void Msg_InputCheck_DofInvalidInElmtBlock(int i,string blockname)
{
    printf("*** Error: %1d-th dof is invalid in [%16s]!      ***\n",i,blockname.c_str());
}