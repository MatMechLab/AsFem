//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"


bool InputSystem::ReadJobBlock(ifstream &in,string str,int &linenum,JobBlock &jobBlock){
    bool HasType=false;
    
    vector<double> numbers;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    str=StrToLower(str);
    numbers.clear();
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            str=StrToLower(str);
            continue;
        }
        str=RemoveStrSpace(str);
        
        if(str.find("type=")!=string::npos||
           str.find("TYPE=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("static")!=string::npos||
               substr.find("STATIC")!=string::npos){
                   jobBlock._JobType=JobType::StaticJob;
                   HasType=true;
            }
            else if(substr.find("transient")!=string::npos||
               substr.find("TRANSIENT")!=string::npos){
                   jobBlock._JobType=JobType::TransientJob;
                   HasType=true;
            }
            else{
                HasType=false;
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported job type in [job] block                  !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        type=static[transient] is expected in [job] block    !!!   ***\n");
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("debug=")!=string::npos||
                str.find("DEBUG=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("true")!=string::npos||
               substr.find("TRUE")!=string::npos){
                jobBlock._IsDebug=true;
                jobBlock._IsDepDebug=false;
            }
            else if(substr.find("false")!=string::npos||
                    substr.find("FALSE")!=string::npos){
                jobBlock._IsDebug=false;
                jobBlock._IsDepDebug=false;
            }
            else if(substr.find("dep")!=string::npos||
                    substr.find("DEP")!=string::npos){
                jobBlock._IsDebug=true;
                jobBlock._IsDepDebug=true;
            }
            else{
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in debug= in [job] block                !!! ***\n");
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("interval=")!=string::npos||
                str.find("INTERVAL=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: interval value not found in [job] block                !!! ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid interval value in [job] block                  !!! ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        interval=integer is expected                           !!! ***\n");
                    Msg_AsFem_Exit();
                }
                jobBlock._Interval=int(numbers[0]);
            }
        }
        else if(str.find("projection=")!=string::npos||
                str.find("PROJECTION=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("true")!=string::npos||
               substr.find("TRUE")!=string::npos){
                jobBlock._IsProjection=true;
            }
            else if(substr.find("false")!=string::npos||
                    substr.find("FALSE")!=string::npos){
                jobBlock._IsProjection=false;
            }
            else{
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in projection= in [job] block           !!! ***\n");
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
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [job] block                          !!! ***\n");
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;
        str=StrToLower(str);
    }
    

    return HasType;
}