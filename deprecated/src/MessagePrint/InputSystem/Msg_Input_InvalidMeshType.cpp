#include "MessagePrint/MessagePrint.h"

void Msg_Input_Invalid1DMeshType()
{
    cout<<"*** Error: invalid 1D mesh type !!!                        ***"<<endl;
    cout<<"***        meshtype=edge2[edge3,edge4] is expected !!!     ***"<<endl;
}

//*******************************+
void Msg_Input_Invalid2DMeshType()
{
    cout<<"*** Error: invalid 2D mesh type !!!                        ***"<<endl;
    cout<<"***        meshtype=quad4,quad8,quad9 is expected !!!      ***"<<endl;
}
//**********************************
void Msg_Input_Invalid3DMeshType()
{
    cout<<"*** Error: invalid 3D mesh type !!!                        ***"<<endl;
    cout<<"***        meshtype=hex8,hex20,hex27 is expected !!!       ***"<<endl;
}