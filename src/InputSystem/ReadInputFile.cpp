//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.30
//+++ Purpose: Function for reading the whole input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadInputFile(Mesh &mesh,DofHandler &dofHandler){
    ifstream in;
    string str;
    int linenum=0;

    bool HasMeshBlock=false;
    bool HasDofsBlock=false;

    if(_HasInputFileName){
        in.open(_InputFileName.c_str(),ios::in);
        while(!in.is_open()){
            MessagePrinter::PrintErrorTxt("can\'t open the input file");
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
            cin>>_InputFileName;
        }
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
            cin>>_InputFileName;
        in.open(_InputFileName.c_str(),ios::in);
        while(!in.is_open()){
            MessagePrinter::PrintErrorTxt("can\'t open the input file");
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the input file name:");
            cin>>_InputFileName;
        }
        _HasInputFileName=true;
    }

    linenum=0;

    HasMeshBlock=false;
    HasDofsBlock=false;
    while(!in.eof()){
        getline(in,str);linenum+=1;
        str=StringUtils::RemoveStrSpace(str);
        str=StringUtils::StrToLower(str);
        if(StringUtils::IsCommentLine(str)||str.size()<1) continue;

        if(str.find("[mesh]")!=string::npos){
            if(!StringUtils::IsBracketMatch(in,linenum)){
                MessagePrinter::PrintErrorTxt("[mesh]/[end] bracket pair dosen\'t match");
                MessagePrinter::AsFem_Exit();
                return false;
            }
            if(ReadMeshBlock(in,str,linenum,mesh)){
                HasMeshBlock=true;
            }
            else{
                HasMeshBlock=false;
            }
        }
        else if(str.find("[dofs]")!=string::npos){
            if(!StringUtils::IsBracketMatch(in,linenum)){
                MessagePrinter::PrintErrorTxt("[dofs]/[end] bracket pair dosen\'t match");
                MessagePrinter::AsFem_Exit();
                return false;
            }
            if(ReadDofsBlock(in,str,linenum,dofHandler)){
                HasDofsBlock=true;
            }
            else{
                HasDofsBlock=false;
            }
        }
    }

    if(!HasMeshBlock){
        MessagePrinter::PrintErrorTxt("no [mesh] block is found, for FEM analysis, mesh is required");
        MessagePrinter::AsFem_Exit();
    }

    if(!HasDofsBlock){
        MessagePrinter::PrintErrorTxt("no [dofs] block is found, for FEM analysis, dofs is required");
        MessagePrinter::AsFem_Exit();
    }

    return true;
}