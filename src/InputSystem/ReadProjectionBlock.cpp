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
//+++ Date   : 2020.07.12
//+++ Purpose: This function can read the [projection] block from our
//+++          input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadProjectionBlock(ifstream &in,string str,int &linenum,SolutionSystem &solutionSystem){
    bool HasName=false;
    vector<string> namelist;
    namelist.clear();
    string str0,str1;
    bool HasScalarName=false;
    bool HasVectorName=false;
    bool HasRank2Name=false;
    bool HasRank4Name=false;
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        getline(in,str);linenum+=1;
        str0=str;
        str1=StringUtils::RemoveStrSpace(str);
        str=StringUtils::StrToLower(str1);
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            continue;
        }
        if(str.compare(0,5,"name=")==0){
            int i=str0.find_first_of('=');
            string substr=str0.substr(i+1,str0.length());
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
        else if(str.compare(0,11,"scalarmate=")==0){
            int i=str0.find_first_of('=');
            string substr=str0.substr(i+1,str0.length());
            namelist=StringUtils::SplitStr(substr,' ');
            if(namelist.size()>0){
                solutionSystem.AddScalarMateProjectionNameFromVec(namelist);
                HasScalarName=true;
            }
            else{
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("scalarmate name can not be found in the [projection] block, 'scalarmate=' should be given");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(str.compare(0,11,"vectormate=")==0){
            int i=str0.find_first_of('=');
            string substr=str0.substr(i+1,str0.length());
            namelist=StringUtils::SplitStr(substr,' ');
            if(namelist.size()>0){
                solutionSystem.AddVectorMateProjectionNameFromVec(namelist);
                HasVectorName=true;
            }
            else{
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("vectormate name can not be found in the [projection] block, 'vectormate=' should be given");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(str.compare(0,10,"rank2mate=")==0){
            int i=str0.find_first_of('=');
            string substr=str0.substr(i+1,str0.length());
            namelist=StringUtils::SplitStr(substr,' ');
            if(namelist.size()>0){
                solutionSystem.AddRank2MateProjectionNameFromVec(namelist);
                HasRank2Name=true;
            }
            else{
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("rank2mate name can not be found in the [projection] block, 'rank2mate=' should be given");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(str.compare(0,10,"rank4mate=")==0){
            int i=str0.find_first_of('=');
            string substr=str0.substr(i+1,str0.length());
            namelist=StringUtils::SplitStr(substr,' ');
            if(namelist.size()>0){
                solutionSystem.AddRank4MateProjectionNameFromVec(namelist);
                HasRank4Name=true;
            }
            else{
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("rank4mate name can not be found in the [projection] block, 'rank4mate=' should be given");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(str.find("[]")!=string::npos&&str.find("[end]")==string::npos){
            MessagePrinter::PrintErrorTxt("the bracket pair is not complete in the [projection] block");
            HasName=false;
            MessagePrinter::AsFem_Exit();
        }
        else if(str.find("[end]")!=string::npos){
            break;
        }
        else{
            MessagePrinter::PrintErrorInLineNumber(linenum);
            MessagePrinter::PrintErrorTxt("unknown option in the [projection] block");
            HasName=false;
            MessagePrinter::AsFem_Exit();
        }
        //getline(in,str);linenum+=1;
    }
    if(!HasName){
        MessagePrinter::PrintErrorTxt("projection name can not be found in the [projection] block, 'name=name1 name2 ...' should be given");
        HasName=false;
        MessagePrinter::AsFem_Exit();
    }

    if(!HasScalarName){
        MessagePrinter::PrintWarningTxt("no scalarmate name found in [projection] block, so the scalar material properties will not be projected");
    }
    if(!HasVectorName){
        MessagePrinter::PrintWarningTxt("no vectormate name found in [projection] block, so the vector material properties will not be projected");
    }
    if(!HasRank2Name){
        MessagePrinter::PrintWarningTxt("no rank2mate name found in [projection] block, so the rank2 tensor material properties will not be projected");
    }
    if(!HasRank4Name){
        MessagePrinter::PrintWarningTxt("no rank4mate name found in [projection] block, so the rank4 tensor material properties will not be projected");
    }

    return HasName;
}