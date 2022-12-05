//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.30
//+++ Purpose: Implement a general class for string manipulate
//+++          For example, string cases convert and the split..
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/StringUtils.h"

StringUtils::StringUtils(){

}

string StringUtils::strToLowerCase(string instr){
    string outstr=instr;
    std::transform(outstr.begin(),outstr.end(),outstr.begin(),
            [](unsigned char c){ return tolower(c);});
    return outstr;
}
string StringUtils::strToUpperCase(string instr){
    string outstr=instr;
    std::transform(outstr.begin(),outstr.end(),outstr.begin(),
            [](unsigned char c){ return toupper(c);});
    return outstr;
}
string StringUtils::removeStrSpace(string instr){
    string outstr=instr;
    outstr.erase(std::remove(outstr.begin(),outstr.end(),' '),outstr.end());
    outstr.erase(std::remove(outstr.begin(),outstr.end(),'\t'),outstr.end());
    outstr.erase(std::remove(outstr.begin(),outstr.end(),'\n'),outstr.end());
    return outstr;
}

vector<double> StringUtils::splitStrNum(string instr){
    vector<double> result;
    int i,n;
    n=instr.size();
    int start,end;
    start=-1;end=-2;
    result.clear();
    if(n==1){
        result.push_back(atof(instr.c_str()));
        return result;
    }
    for(i=0;i<n;i++){
        if(i==0){
            if(instr.at(i)=='-'&&(instr.at(i+1)>='0'&&instr.at(i+1)<='9')){
                start = i;
            }
            else if(instr.at(i)>='0'&&instr.at(i)<='9'){
                start=i;
            }
            if((instr.at(i)>='0'&&instr.at(i)<='9')&&(instr.at(i+1)==' '||instr.at(i+1)==',')){
                end=i;
            }
        }
        else if(i==n-1){
            if(instr.at(i)>='0'&&instr.at(i)<='9'){
                end=i;
                if((instr.at(i-1)<'0'||instr.at(i-1)>'9')&&instr.at(i-1)!='.'&&instr.at(i-1)!='e'&&instr.at(i-1)!='-'){
                    start=i;
                }
            }
        }
        else{
            if(instr.at(i)=='-'&&(instr.at(i+1)>='0'&&instr.at(i+1)<='9')){
                if(instr.at(i-1)==' '||instr.at(i-1)==','||instr.at(i-1)=='='||instr.at(i-1)=='<'||instr.at(i-1)=='>'||instr.at(i-1)==':'){
                    start=i;
                }
            }
            else if((instr.at(i)>='0'&&instr.at(i)<='9')&&((instr.at(i-1)<'0'||instr.at(i-1)>'9')&&instr.at(i-1)!='.'&&instr.at(i-1)!='e'&&instr.at(i-1)!='-')){
                start=i;
            }

            if((instr.at(i)>='0'&&instr.at(i)<='9')&&(instr.at(i+1)<'0'||instr.at(i+1)>'9')&&instr.at(i+1)!='.'&&instr.at(i+1)!='e'){
                end=i;
            }
        }
        if(end>=start&&end==i&&start>=0){
            result.push_back(atof(instr.substr(start,end+1-start).c_str()));
        }
    }
    return result;
}

bool StringUtils::isValidTimeDependentExpression(string instr){
    // now the string should contains:
    // "1.0*t" or "-1.0*t"
    string str=removeStrSpace(instr);
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
    if(splitStrNum(str1).size()>0&&str2.compare(0,2,"*t")==0){
        return true;
    }
    // for t*number case
    str1=str.substr(0,i+1);
    str2=str.substr(i+1);
    if(str1.compare(0,2,"t*")==0&&splitStrNum(str2).size()>0){
        return true;
    }
    return false;
}