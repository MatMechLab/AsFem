#include "MessagePrint/MessagePrint.h"

void Msg_Input_LinearSolverBlock_InvalidSolverOption(){
    cout<<"*** Error: invalid solver option in [linearsolver] block !!***"<<endl;
    cout<<"***        solver=lu[cg,bicg,qr...] is expected !!!        ***"<<endl;
}