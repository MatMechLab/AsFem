#include "MessagePrint/MessagePrint.h"

void Msg_GaussPoint_Invalid3DMeshType()
{
    printf("*** Error: unsupported mesh type for 3D gauss point!!!     ***\n");
    printf("***        hex8,20,27/tet4,10 is expected for 3D case!     ***\n");
}