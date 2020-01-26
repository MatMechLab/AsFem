#include "MessagePrint/MessagePrint.h"

void Msg_Input_LinearSolverBlock_TolNotFound(){
    cout<<"*** Error: tolerance not found in [linearsolver] block !!! ***"<<endl;
    cout<<"***        tol=some positive real value is expected !!!    ***"<<endl;
}