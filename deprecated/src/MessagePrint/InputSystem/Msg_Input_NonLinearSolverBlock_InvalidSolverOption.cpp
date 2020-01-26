#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_InvalidSolverOption(){
    cout<<"*** Error: invalid type option in [nonlinearsolver] block !***"<<endl;
    cout<<"***        type=nr[nrmodify,nrlinesearch,...] is expected !***"<<endl;
}