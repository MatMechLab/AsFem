//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadDofsBlock(ifstream &in,string str,int &linenum,DofHandler &dofHandler){
    // dof block format:
    // [dofs]
    //   name=c1 c2
    // [end]
    bool HasName=false;
    vector<string> namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    // str=StrToLower(str);
    namelist.clear();
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            // str=StrToLower(str);
            continue;
        }
        if(str.find("name=")!=string::npos||
           str.find("NAME=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            namelist=SplitStr(substr,' ');
            if(namelist.size()>0)
            {
                if(IsUniqueStrVec(namelist)){
                    HasName=true;
                    dofHandler.AddNameToDofNameList(namelist);
                }
                else{
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: duplicated dof name in [dofs] block !!!                    ***\n");
                    Msg_AsFem_Exit();
                }
            }
            else{
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dof found in [dofs] block                !!!            ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        name=dof1 dof2 ... is expected              !!!            ***\n");
                HasName=false;
                Msg_AsFem_Exit();
            }
            
        }
        else if(str.find("[]")!=string::npos){
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
        else{
            Msg_Input_LineError(linenum);
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [dofs] block              !!!            ***\n");
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;
        str=StrToLower(str);
    }
    return HasName;
}