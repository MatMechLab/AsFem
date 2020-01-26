#include "MessagePrint/MessagePrint.h"

void Msg_GaussPoint_Invalid3DLobattoMeshType()
{
    printf("*** Error: unsupported mesh type for 3D gauss point!!!     ***\n");
    printf("***        hex8,20,27 is expected for 3D lobatto case!     ***\n");
}