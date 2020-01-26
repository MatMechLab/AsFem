#include "MessagePrint/MessagePrint.h"

void Msg_Input_JobBlock_InvalidSolverOption(){
    cout<<"*** Error: invalid solver option in your input file !!     ***"<<endl;
    cout<<"***        solver=lu[cg,bicg,qr...] is expected !!!        ***"<<endl;
}