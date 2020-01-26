#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_ErelTolInvalid(){
    cout<<"*** Error: invalid e_rel_tol in [nonlinearsolver] block !!!***"<<endl;
    cout<<"***        e_rel_tol=positive real value is expected !!!   ***"<<endl;
}