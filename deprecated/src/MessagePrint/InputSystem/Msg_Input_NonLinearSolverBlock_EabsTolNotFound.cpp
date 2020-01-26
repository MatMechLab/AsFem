#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_EabsTolNotFound(){
    cout<<"*** Error: e_abs_tol not found in [nonlinearsolver] block !***"<<endl;
    cout<<"***        e_abs_tol=positive real value is expected !!!   ***"<<endl;
}