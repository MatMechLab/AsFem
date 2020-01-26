#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_EabsTolInvalid(){
    cout<<"*** Error: invalid e_abs_tol in [nonlinearsolver] block !!!***"<<endl;
    cout<<"***        e_abs_tol=positive real value is expected !!!   ***"<<endl;
}