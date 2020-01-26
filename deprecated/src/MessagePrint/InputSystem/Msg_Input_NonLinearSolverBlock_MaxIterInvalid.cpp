#include "MessagePrint/MessagePrint.h"

void Msg_Input_NonLinearSolverBlock_MaxIterInvalid(){
    cout<<"*** Error: invalid maxiters in [nonlinearsolver] block !!! ***"<<endl;
    cout<<"***        maxiters=some positive integer is expected !!!  ***"<<endl;
}