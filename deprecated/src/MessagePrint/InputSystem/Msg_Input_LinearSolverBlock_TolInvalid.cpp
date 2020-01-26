#include "MessagePrint/MessagePrint.h"

void Msg_Input_LinearSolverBlock_TolInvalid(){
    cout<<"*** Error: invalid tolerance in [linearsolver] block !!!   ***"<<endl;
    cout<<"***        tol=some positive real value is expected !!!    ***"<<endl;
}