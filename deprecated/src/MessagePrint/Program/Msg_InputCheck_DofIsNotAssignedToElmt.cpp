#include "MessagePrint/MessagePrint.h"

void Msg_InputCheck_DofIsNotAssignedToElmt(string dofname){
    printf("*** Error: dof=%-10s isn\'t assigned to any elmts!     ***\n",dofname.c_str());
    printf("***        you should apply an elmt to it !!!              ***\n");
}