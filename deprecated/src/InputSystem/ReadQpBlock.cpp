#include "InputSystem/InputSystem.h"

bool InputSystem::ReadQpBlock(ifstream &in,string str,int &linenum,QPBlockInfo &qpBlockInfo)
{
    // dof block format:
    // [qpoint]
    //   type=gauss [gausslobatto]
    //   order=4
    // [end]
    

    bool HasType=false;
    bool HasOrder=false;
    vector<double> numbers;
    string namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;

    qpBlockInfo._nQpOrder=1;// if no order is given, then default order is 1
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos)
    {
        if(IsCommentLine(str)||str.length()<1)
        {
            getline(in,str);linenum+=1;
            continue;
        }
        if(str.find("type=")!=string::npos||
           str.find("TYPE=")!=string::npos)
        {
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if((substr.find("gauss")!=string::npos||substr.find("GAUSS")!=string::npos)&&
               substr.length()==5){
                HasType=true;
                qpBlockInfo._SetFromInput=true;
                qpBlockInfo._QpType="gauss";
            }
            else if((substr.find("gausslobatto")!=string::npos||
                     substr.find("GAUSSLOBATTO")!=string::npos)&&
                    substr.length()==12){
                HasType=true;
                qpBlockInfo._SetFromInput=true;
                qpBlockInfo._QpType="gausslobatto";
            }
            else
            {
                Msg_Input_LineError(linenum);
                Msg_Input_QpBlockTypeNotFound();
                HasType=false;
                Msg_AsFem_Exit();
            }
            
        }
        else  if(str.find("order")!=string::npos||str.find("ORDER")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [qpoint] block  !!!        ***"<<endl;
                cout<<"***        order should be given after type in [qpoint] !!!***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_QpBlockOrderNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<0||int(numbers[0])>11){
                    Msg_Input_LineError(linenum);
                    Msg_Input_QpBlockOrderInvalid();
                    Msg_AsFem_Exit();
                }
                qpBlockInfo._nQpOrder=int(numbers[0]);
                HasOrder=true;
            }
        }
        else if(str.find("[]")!=string::npos)
        {
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
        else{
            Msg_Input_LineError(linenum);
            cout<<"*** Error: unknown option in [projection] block    !!!     ***"<<endl;
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;
    }
    if(!HasOrder){
        cout<<"*** Warning: no order is asigned in [qpoint] block !!!     ***"<<endl;
        cout<<"***          order=1 will be used for integration !!!      ***"<<endl;
    }

    return HasType;
}