#include "MessagePrint/MessagePrint.h"

void Msg_GaussPoint_Invalid2DLobattoMeshType()
{
    printf("*** Error: unsupported mesh type for 2D gauss point!!!     ***\n");
    printf("***        quad4,8,9 is expected for 2D lobatto case!!     ***\n");
}