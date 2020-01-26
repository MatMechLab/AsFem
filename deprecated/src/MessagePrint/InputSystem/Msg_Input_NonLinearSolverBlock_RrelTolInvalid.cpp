#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_RrelTolInvalid(){
    cout<<"*** Error: invalid r_rel_tol in [nonlinearsolver] block!!! ***"<<endl;
    cout<<"***        r_rel_tol=positive real value is expected !!!   ***"<<endl;
}