#include "MessagePrint/MessagePrint.h"

void Msg_GaussPoint_Invalid3DGaussOrder(int n)
{
    printf("*** Error: unsupported gauss integration order(=%2d)!!     ***\n",n);
    printf("***        n=0~3 is expected for 2D case !!!               ***\n");
}