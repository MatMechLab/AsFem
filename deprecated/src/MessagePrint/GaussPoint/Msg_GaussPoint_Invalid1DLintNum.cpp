#include "MessagePrint/MessagePrint.h"

void Msg_GaussPoint_Invalid1DLintNum(int n)
{
    printf("*** Error: unsupported gauss integration num(=%2d) !!!     ***\n",n);
    printf("***        n=1~5 is expected for 1D case !!!               ***\n");
}