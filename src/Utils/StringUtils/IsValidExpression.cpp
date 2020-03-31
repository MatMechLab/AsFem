//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/StringUtils.h"

bool IsValidExpression(string str){
    // now the string should contains:
    // "1.0*t" or "-1.0*t"
    str=RemoveStrSpace(str);
    if(str.size()<1){
        return false;
    }
    else if(str.size()==1){
        if(str.find("t")!=string::npos){
            return true;
        }
        else{
            return false;
        }
    }
    else if(str.size()==2){
        if(str.find("-t")!=string::npos){
            return true;
        }
        else{
            return false;
        }
    }

    if(str.find("t*")==string::npos&&str.find("*t")==string::npos){
        return false;
    }
    int i=str.find_first_of('*');

    // for number*t case
    string str1=str.substr(0,i);
    string str2=str.substr(i);
    // cout<<"str1="<<str1<<", str2="<<str2<<endl;
    if(SplitStrNum(str1).size()>0&&str2.compare(0,2,"*t")==0){
        return true;
    }
    // for t*number case
    str1=str.substr(0,i+1);
    str2=str.substr(i+1);
    // cout<<"str1="<<str1<<", str2="<<str2<<endl;
    if(str1.compare(0,2,"t*")==0&&SplitStrNum(str2).size()>0){
        return true;
    }
    return false;
}