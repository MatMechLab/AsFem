#include "MessagePrint/MessagePrint.h"

void Msg_GaussPoint_Invalid1DLobbotLintNum(int n)
{
    printf("*** Error: unsupported gauss integration num(=%2d) !!!     ***\n",n);
    printf("***        n=3~6 is expected for 1D lobbot case !!!        ***\n");
}