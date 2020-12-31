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
//+++ Date   : 2020.07.01
//+++ Purpose: This function can read the [dofs] block from our
//+++          input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadDofsBlock(ifstream &in,string str,int &linenum,DofHandler &dofHandler){
    // dof block format:
    // [dofs]
    //   name=c1 c2
    // [end]
    char buff[55];
    bool HasName=false;
    vector<string> namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    // str=StrToLower(str);
    namelist.clear();
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            // str=StrToLower(str);
            continue;
        }
        if(str.find("name=")!=string::npos||
           str.find("NAME=")!=string::npos){
            int i=str.find_first_of('=');
            if(str.size()<=5){
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt(" no dof name found in the [dofs] block, name=dof1 dof2 ... is expected");
                MessagePrinter::AsFem_Exit();
            }
            string substr=str.substr(i+1,str.length());
            namelist=StringUtils::SplitStr(substr,' ');
            if(namelist.size()>0){
                if(StringUtils::IsUniqueStrVec(namelist)){
                    HasName=true;
                    dofHandler.AddDofNameFromStrVec(namelist);
                }
                else{
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt(" duplicated dof name in the [dofs] block");
                    MessagePrinter::AsFem_Exit();
                }
            }
            else{
                HasName=false;

                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt(" no dof name found in the [dofs] block, name=dof1 dof2 ... is expected");
                MessagePrinter::AsFem_Exit();
            }
            
        }
        else if(str.find("[]")!=string::npos){
            snprintf(buff,55,"line-%d has some errors",linenum);
            MessagePrinter::PrintErrorTxt(string(buff));
            MessagePrinter::PrintErrorTxt(" the bracket pair is not complete in the [dofs] block");
            MessagePrinter::AsFem_Exit();
        }
        else{
            snprintf(buff,55,"line-%d has some errors",linenum);
            MessagePrinter::PrintErrorTxt(string(buff));
            MessagePrinter::PrintErrorTxt(" unknown option in the [dofs] block");
            MessagePrinter::AsFem_Exit();
        }
        getline(in,str);linenum+=1;
        str=StringUtils::StrToLower(str);
    }
    if(!HasName){
        MessagePrinter::PrintErrorTxt(" no dof name found in the [dofs] block, name=dof1 dof2 ... is expected");
        MessagePrinter::AsFem_Exit();
    }
    return HasName;
}