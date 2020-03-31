//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadQPointBlock(ifstream &in,string str,int &linenum,FE &fe){
    // format:
    // [qpoint]
    //   type=gauss[gausslobatto]
    //   order=2
    //   bcorder=1
    // [end]

    bool HasType;
    vector<double> numbers;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    str=StrToLower(str);
    numbers.clear();

    HasType=false;
    
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
            if((substr.find("gauss")!=string::npos||
               substr.find("GAUSS")!=string::npos)&&substr.size()==5){
                   fe.SetQPointType("gauss");
                   HasType=true;
            }
            else if((substr.find("gausslobatto")!=string::npos||
                    substr.find("GAUSSLOBATTO")!=string::npos)&&substr.size()==12){
                   fe.SetQPointType("gausslobatto");
                   HasType=true;
            }
            else{
                HasType=false;
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported gauss point type in [qpoint] block       !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        type=gauss[gausslobatto] is expected in [qpoint]block!!!   ***\n");
                Msg_AsFem_Exit();
            }
        }
        else if((str.find("order=")!=string::npos||
                 str.find("ORDER=")!=string::npos)&&
                 str.find("bcorder")==string::npos&&
                 str.find("BCORDER")==string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: order= not found in [qpoint] block                     !!! ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<0||int(numbers[0])>7){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid order value in [qpoint] block                  !!! ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        order=integer is expected                              !!! ***\n");
                    Msg_AsFem_Exit();
                }
                fe.SetOrder(int(numbers[0]));
                fe.SetBCOrder(int(numbers[0]));
            }
        }
        else if(str.find("bcorder=")!=string::npos||
                str.find("BCORDER=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: bcorder= not found in [qpoint] block                   !!! ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<0||int(numbers[0])>7){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid bcorder value in [qpoint] block                !!! ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        bcorder=integer is expected                            !!! ***\n");
                    Msg_AsFem_Exit();
                }
                fe.SetBCOrder(int(numbers[0]));
            }
        }
        else if(str.find("[]")!=string::npos){
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
        else{
            Msg_Input_LineError(linenum);
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [qpoint] block                       !!! ***\n");
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;
        str=StrToLower(str);
    }
    

    return HasType;
}