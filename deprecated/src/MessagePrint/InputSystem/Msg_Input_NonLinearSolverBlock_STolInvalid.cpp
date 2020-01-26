#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_STolInvalid(){
    cout<<"*** Error: invalid s_tol in [nonlinearsolver] block !!!    ***"<<endl;
    cout<<"***        s_tol=positive real value is expected !!!       ***"<<endl;
}