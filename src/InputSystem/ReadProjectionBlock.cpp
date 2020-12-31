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
//+++ Date   : 2020.07.12
//+++ Purpose: This function can read the [projection] block from our
//+++          input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadProjectionBlock(ifstream &in,string str,int &linenum,SolutionSystem &solutionSystem){
    bool HasName=false;
    vector<string> namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    namelist.clear();
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            continue;
        }
        if(str.find("name=")!=string::npos||
           str.find("NAME=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            namelist=StringUtils::SplitStr(substr,' ');
            if(namelist.size()>0){
                solutionSystem.AddProjectionNameFromVec(namelist);
                HasName=true;
            }
            else{
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("projection name can not be found in the [projection] block");
                HasName=false;
                MessagePrinter::AsFem_Exit();
            }
            
        }
        else if(str.find("[]")!=string::npos){
            MessagePrinter::PrintErrorTxt("the bracket pair is not complete in the [projection] block");
            HasName=false;
            MessagePrinter::AsFem_Exit();
        }
        else{
            MessagePrinter::PrintErrorInLineNumber(linenum);
            MessagePrinter::PrintErrorTxt("unknown option in the [projection] block");
            HasName=false;
            MessagePrinter::AsFem_Exit();
        }
        getline(in,str);linenum+=1;
    }
    if(!HasName){
        MessagePrinter::PrintErrorTxt("projection name can not be found in the [projection] block, 'name=name1 name2 ...' should be given");
        HasName=false;
        MessagePrinter::AsFem_Exit();
    }

    return HasName;
}