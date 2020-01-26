#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_RrelTolNotFound(){
    cout<<"*** Error: r_rel_tol not found in [nonlinearsolver] block!!***"<<endl;
    cout<<"***        r_rel_tol=positive real value is expected !!!   ***"<<endl;
}