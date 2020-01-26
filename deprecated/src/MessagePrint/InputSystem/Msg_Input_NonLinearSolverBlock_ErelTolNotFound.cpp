#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_ErelTolNotFound(){
    cout<<"*** Error: e_rel_tol not found in [nonlinearsolver] block !***"<<endl;
    cout<<"***        e_rel_tol=positive real value is expected !!!   ***"<<endl;
}