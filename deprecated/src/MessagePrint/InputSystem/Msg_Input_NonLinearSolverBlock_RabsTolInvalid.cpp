#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_RabsTolInvalid(){
    cout<<"*** Error: invalid r_abs_tol [nonlinearsolver] block !!!   ***"<<endl;
    cout<<"***        r_abs_tol=positive real value is expected !!!   ***"<<endl;
}