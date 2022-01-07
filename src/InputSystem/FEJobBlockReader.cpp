//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.08.06
//+++ Purpose: Implement the reader for [job] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "InputSystem/FEJobBlockReader.h"

void FEJobBlockReader::PrintHelper(){
    MessagePrinter::PrintStars(MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("The complete information for [job] block:",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[job]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  type=static,transient",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  debug=true,false,dep",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[end]",MessageColor::BLUE);
    MessagePrinter::PrintStars(MessageColor::BLUE);
}
//****************************************************************
bool FEJobBlockReader::ReadFEJobBlock(ifstream &in, string str, int &linenum, FEJobBlock &feJobBlock){
    char buff[55];
    bool HasType=false;
    vector<string> namelist;
    // now str already contains [job]
    getline(in,str);linenum+=1;
    str=StringUtils::RemoveStrSpace(str);
    str=StringUtils::StrToLower(str);

    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            str=StringUtils::RemoveStrSpace(str);
            str=StringUtils::StrToLower(str);
            continue;
        }
        if(str.find("type=")!=string::npos||
           str.find("TYPE=")!=string::npos){
            int i=str.find_first_of('=');
            if(str.size()<=5){
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt(" no type name found in the [job] block, type=static[transient] is expected");
                MessagePrinter::AsFem_Exit();
            }
            string substr=str.substr(i+1,str.length());
            if(substr.find("helper")!=string::npos){
                PrintHelper();
                return false;
            }
            else if(substr.find("static")!=string::npos||
               substr.find("STATIC")!=string::npos){
                feJobBlock._jobType=FEJobType::STATIC;
                feJobBlock._jobTypeName="static";
                HasType=true;
            }
            else if(substr.find("transient")!=string::npos||
               substr.find("TRANSIENT")!=string::npos){
                feJobBlock._jobType=FEJobType::TRANSIENT;
                feJobBlock._jobTypeName="transient";
                HasType=true;
            }
            else{
                HasType=false;
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt("unsupported type name found in the [job] block, type=static[transient] is expected");
                MessagePrinter::AsFem_Exit();
            }

        }
        else if(str.find("debug=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            if(substr.find("false")!=string::npos){
                feJobBlock._IsDebug=false;
                feJobBlock._IsDepDebug=false;
            }
            else if(substr.find("true")!=string::npos||
               substr.find("TRUE")!=string::npos){
                feJobBlock._IsDebug=true;
                feJobBlock._IsDepDebug=false;
            }
            else if(substr.find("dep")!=string::npos||
                substr.find("DEP")!=string::npos){
                feJobBlock._IsDebug=true;
                feJobBlock._IsDepDebug=true;
            }
            else{
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt(" unknown debug option for debug= in [job] block");
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
        str=StringUtils::RemoveStrSpace(str);
    }
    if(!HasType){
        MessagePrinter::PrintErrorTxt(" no type name found in the [job] block, type=static[transient] is expected");
        MessagePrinter::AsFem_Exit();
    }
    return HasType;
}
