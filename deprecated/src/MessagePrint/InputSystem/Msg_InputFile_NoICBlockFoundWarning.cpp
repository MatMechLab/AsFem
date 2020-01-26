#include "MessagePrint/MessagePrint.h"

void Msg_InputFile_NoICBlockFoundWarning()
{
    cout<<"*** Warning: no any initial conditions are given !!!       ***"<<endl;
    cout<<"***          all the dofs will be set to zero !!!          ***"<<endl;
}