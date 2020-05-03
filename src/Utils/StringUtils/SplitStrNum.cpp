//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/StringUtils.h"

vector<double> SplitStrNum(string &instr)
{
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
        //cout<<"i="<<i<<" start="<<start<<"<--->end="<<end<<endl;
        if(end>=start&&end==i&&start>=0){
            //cout<<"start="<<start<<"<--->end="<<end<<endl;
            result.push_back(atof(instr.substr(start,end+1-start).c_str()));
        }
    }
    return result;
}

vector<double> SplitStrNum(string &instr,char symbol){
    // cout<<"instr="<<instr<<endl;
    string tempstr;
    tempstr=instr;
    if(symbol!=' '){
        tempstr=RemoveStrSpace(instr);
    }
    else{
        tempstr=instr;
    }
    // cout<<"str size="<<tempstr.size()<<", str length="<<tempstr.length()<<endl;
    // for(int i=0;i<static_cast<int>(tempstr.length());i++){
    //     cout<<"i="<<i+1<<":"<<tempstr.at(i)<<endl;
    // }
    vector<string> strlist=SplitStr(tempstr,symbol);
    // for(int i=0;i<static_cast<int>(strlist.size());i++){
    //     cout<<"i="<<i+1<<", size="<<strlist[i].size()<<", content="<<strlist[i]<<endl;
    // }
    if(strlist.size()>1){
        if(strlist[strlist.size()-1].size()==1){
            if(strlist[strlist.size()-1].at(0)<'0'||strlist[strlist.size()-1].at(0)>'9'){
                strlist.pop_back();
            }
        }
    }
    // for(auto it:strlist)cout<<it<<endl;
    vector<double> numbers;
    numbers.clear();
    for(int i=0;i<static_cast<int>(strlist.size());i++){
        // cout<<"strchar="<<strlist[i].c_str()<<endl;
        numbers.push_back(stod(strlist[i]));
    }
    return numbers;
}