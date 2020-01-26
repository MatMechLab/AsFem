#include "MessagePrint/MessagePrint.h"

void Msg_Input_JobBlock_InvalidNonLinearSolverOption(){
    cout<<"*** Error: invalid nonlinear method in your input file!    ***"<<endl;
    cout<<"***        method=nr[modifynr,...] is expected !!!         ***"<<endl;
}