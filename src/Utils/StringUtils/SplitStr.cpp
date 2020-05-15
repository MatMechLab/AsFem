//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/StringUtils.h"

vector<string> SplitStr(string &instr,char symbol)
{
    string instr0=instr;
    vector<string> outstr;
    outstr.clear();

    if(instr0.find_last_of('\n')!=string::npos){
        instr0.pop_back();
    }
    else if(instr0.at(instr0.size()-1)==symbol){
        instr0.pop_back();
    }

    stringstream ss(instr0);
    string temp;
    while (getline(ss,temp,symbol)){
        outstr.push_back(temp);
    }
    vector<string> outstr0;
    outstr0.clear();
    if(symbol==' '){
        for(unsigned int i=0;i<outstr.size();i++){
            temp=outstr[i];
            temp.erase(remove(temp.begin(),temp.end(),' '),temp.end());
            if(temp.length()>0){
                outstr0.push_back(temp);
            }
        }
        outstr=outstr0;
    }
    else{
        outstr0=outstr;
    }
    outstr=outstr0;
    for(auto &it:outstr){
        it.erase(remove(it.begin(),it.end(),symbol),it.end());
        it.erase(remove(it.begin(),it.end(),'\n'),it.end());
    }
    return outstr;
}