//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadOutputBlock(ifstream &in,string str,int &linenum,OutputBlock &outputblock){
    // format:
    // [output]
    //   type=vtu[vtk,csv,txt]
    //   folder=foldername[default is empty]
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
            if((substr.find("vtu")!=string::npos||
               substr.find("VTU")!=string::npos)&&substr.size()==3){
                   outputblock._OutputFileType=OutputFileType::VTU;
                   HasType=true;
            }
            // else if((substr.find("vtk")!=string::npos||
            //    substr.find("VTK")!=string::npos)&&substr.size()==3){
            //        outputblock._OutputFileType=OutputFileType::VTK;
            //        HasType=true;
            // }
            // else if((substr.find("csv")!=string::npos||
            //    substr.find("CSV")!=string::npos)&&substr.size()==3){
            //        outputblock._OutputFileType=OutputFileType::CSV
            //        HasType=true;
            // }
            else{
                HasType=false;
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported output file type in [output] block       !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        type=vtu[vtk,csv] is expected in [output] block      !!!   ***\n");
                Msg_AsFem_Exit();
            }
        }
        else if((str.find("folder=")!=string::npos||
                 str.find("Folder=")!=string::npos)&&
                 str.find("FOLDER")==string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.size()<1){
                outputblock._FolderName.clear();
            }
            else{
                outputblock._FolderName=substr;
            }
        }
        else if(str.find("[]")!=string::npos){
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
        else{
            Msg_Input_LineError(linenum);
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [output] block                       !!! ***\n");
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;
        str=StrToLower(str);
    }

    HasType=true;
    

    return HasType;
}