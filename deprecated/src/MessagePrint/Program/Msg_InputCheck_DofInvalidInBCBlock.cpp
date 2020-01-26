#include "MessagePrint/MessagePrint.h"


void Msg_InputCheck_DofInvalidInBCBlock(string blockname)
{
    printf("*** Error: dof name is invalid in [%16s] !     ***\n",blockname.c_str());
}