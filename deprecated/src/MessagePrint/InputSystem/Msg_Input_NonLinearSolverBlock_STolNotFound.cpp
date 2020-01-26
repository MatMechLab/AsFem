#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_STolNotFound(){
    cout<<"*** Error: s_tol not found in [nonlinearsolver] block!!    ***"<<endl;
    cout<<"***        s_tol=positive real value is expected !!!       ***"<<endl;
}