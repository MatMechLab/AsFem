#include "MessagePrint/MessagePrint.h"

void Msg_GaussPoint_Invalid2DMeshType()
{
    printf("*** Error: unsupported mesh type for 2D gauss point!!!     ***\n");
    printf("***        quad4,8,9/tri3,6 is expected for 2D case!!!     ***\n");
}