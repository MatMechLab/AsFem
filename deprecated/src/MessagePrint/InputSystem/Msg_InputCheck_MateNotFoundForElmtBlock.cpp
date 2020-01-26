#include "MessagePrint/MessagePrint.h"

void Msg_InputCheck_MateNotFoundForElmtBlock(string matename,string blockname)
{
    if(matename.size()<1){
        printf("*** Error: no mate name found for [%13s] !!!      ***\n",blockname.c_str());
        cout<<"***        [mates] block is required for [elmts] block!!!  ***"<<endl;
    }
    else{
        printf("*** Error: no [%11s] found for [%13s]!     ***\n",matename.c_str(),blockname.c_str());
        printf("***        [%11s] is required in [mates] block!     ***\n",matename.c_str());
    }
}