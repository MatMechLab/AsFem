#include "InputSystem/InputSystem.h"

bool InputSystem::ReadProjectionBlock(ifstream &in,string str,int &linenum,SolutionSystem &solutionSystem){
    // projection block format:
    // [projection]
    //   name=v1 v2 v3 ...
    // [end
    

    bool HasName=false;
    vector<string> namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    namelist.clear();
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            continue;
        }
        if(str.find("name=")!=string::npos||
           str.find("NAME=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            namelist=SplitStr(substr,' ');
            if(namelist.size()>0){
                solutionSystem.SetProjNameFromStrVec(namelist);
                HasName=true;
            }
            else{
                Msg_Input_LineError(linenum);
                Msg_Input_ProjectionNameNotFound();
                HasName=false;
                Msg_AsFem_Exit();
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
    if(!HasName){
        cout<<"*** Error: no 'name=' found in [projection] block !!!      ***"<<endl;
        cout<<"***        name=x y z should be given !!!                  ***"<<endl;
        Msg_AsFem_Exit();
    }

    return HasName;
}