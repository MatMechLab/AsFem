#include "InputSystem/InputSystem.h"

bool InputSystem::ReadDofsBlock(ifstream &in,string str,int &linenum,DofHandler &dofHandler)
{
    // dof block format:
    // [dofs]
    //   name=c1 c2
    // [end]
    

    bool HasName=false;
    vector<string> namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    namelist.clear();
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos)
    {
        if(IsCommentLine(str)||str.length()<1)
        {
            getline(in,str);linenum+=1;
            continue;
        }
        if(str.find("name=")!=string::npos||
           str.find("NAME=")!=string::npos)
        {
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            namelist=SplitStr(substr,' ');
            if(namelist.size()>0)
            {
                if(IsUniqueStrVec(namelist)){
                    HasName=true;
                    dofHandler.SetDofNameListFromVec(namelist);
                }
                else{
                    Msg_Input_LineError(linenum);
                    Msg_Input_DofNameDuplicate();
                    Msg_AsFem_Exit();
                }
            }
            else
            {
                Msg_Input_LineError(linenum);
                Msg_Input_DofNameNotFound();
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
        getline(in,str);linenum+=1;
    }

    return HasName;
}