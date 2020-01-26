#include "MessagePrint/MessagePrint.h"

void Msg_Input_JobBblock_TransientInfoNotComplete(){
    cout<<"*** Error: transient job information not complete!!!       ***"<<endl;
    cout<<"***        endtime and dt must be given in [job]!!!        ***"<<endl;
}