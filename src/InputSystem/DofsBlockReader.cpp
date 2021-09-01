//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.07.13
//+++ Purpose: Implement the reader for [dofs] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/DofsBlockReader.h"

void DofsBlockReader::PrintHelper(){
    MessagePrinter::PrintStars(MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("The complete information for [dofs] block:",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[dofs]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("name=dof1_name dof2_name ...",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[end]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("each name should be separated by a space",MessageColor::BLUE);
    MessagePrinter::PrintStars(MessageColor::BLUE);
}

//*************************************************************
bool DofsBlockReader::ReadDofsBlock(ifstream &in, string str, int &linenum, DofHandler &dofHandler){
    char buff[55]; /**< this char vector used to convert message to string */
    bool HasName=false;
    vector<string> namelist; /**< store all the single dof name from input */
    string str0;
    // now the str already contains '[dofs]'
    getline(in,str);linenum+=1;
    namelist.clear();

    str0=str;
    str=StringUtils::StrToLower(str0);
    str=StringUtils::RemoveStrSpace(str);
    
    while(str.find("[end]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            str0=str;
            str=StringUtils::StrToLower(str);
            str=StringUtils::RemoveStrSpace(str);
            continue;
        }
        if(str.find("type=helper")!=string::npos){
            PrintHelper();
            return false;
        }
        else if(str.find("name=")!=string::npos){
            int i=str0.find_first_of('=');
            if(str.size()<=5){
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt(" no dof name found in the [dofs] block, name=dof1 dof2 ... is expected");
                MessagePrinter::AsFem_Exit();
            }

            string substr=str0.substr(i+1,str0.length());
            namelist=StringUtils::SplitStr(substr,' ');/**< here we use the original string to keep the upper and low case for Dofs name */
            if(namelist.size()>0){
                if(StringUtils::IsUniqueStrVec(namelist)){
                    /**
                     * we only add the unique name to our dofHandler class, otherwise
                     * the input file is invalid!
                     */
                    HasName=true;
                    dofHandler.AddDofNameFromStrVec(namelist);
                    //for(auto it:namelist) cout<<it<<endl;
                }
                else{
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt(" duplicated dof name in the [dofs] block");
                    MessagePrinter::AsFem_Exit();
                } /**< check either the name is unique or not */
            }
            else{
                /**< if no dof name is found ! */
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
            /**< for some unknown options in [dofs] block */
            snprintf(buff,55,"line-%d has some errors",linenum);
            MessagePrinter::PrintErrorTxt(string(buff));
            MessagePrinter::PrintErrorTxt(" unknown option in the [dofs] block");
            MessagePrinter::AsFem_Exit();
        }
        getline(in,str);linenum+=1;
    
    } /**< end of while loop */
    
    return HasName;
}
