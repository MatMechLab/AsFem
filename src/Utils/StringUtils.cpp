//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
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

string StringUtils::StrToLower(string instr){
    string outstr=instr;
    transform(outstr.begin(),outstr.end(),outstr.begin(),
    [](unsigned char c){ return tolower(c);});
    return outstr;
}
string StringUtils::StrToUpper(string instr){
    string outstr=instr;
    transform(outstr.begin(),outstr.end(),outstr.begin(),
    [](unsigned char c){ return toupper(c);});
    return outstr;
}
//************************************************
string StringUtils::RemoveStrSpace(string &instr){
    string outstr=instr;
    outstr.erase(remove(outstr.begin(),outstr.end(),' '),outstr.end());
    outstr.erase(remove(outstr.begin(),outstr.end(),'\t'),outstr.end());
    outstr.erase(remove(outstr.begin(),outstr.end(),'\n'),outstr.end());
    return outstr;
}
string StringUtils::RemoveSymbolFromStr(string &instr,char symbol){
    string outstr=instr;
    outstr.erase(remove(outstr.begin(),outstr.end(),symbol),outstr.end());
    return outstr;
}
vector<string> StringUtils::SplitStr(string &instr,char symbol){
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
        it.erase(remove(it.begin(),it.end(),'\t'),it.end());
    }
    return outstr;
}
//********************************************************
bool StringUtils::IsUniqueStrVec(vector<string> &strvec){
    if(strvec.size()<=1){
        return true;
    }
    else if(strvec.size()==2){
        if(strvec[0]==strvec[1]){
            return false;
        }
        else{
            return true;
        }
    }
    else if(strvec.size()==3){
        if(strvec[0]!=strvec[1]&&strvec[0]!=strvec[2]&&strvec[1]!=strvec[2]){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        bool IsDuplicate=false;
        for(unsigned int i=0;i<strvec.size();++i){
            for(unsigned int j=0;j<strvec.size();++j){
                if(i!=j){
                    if(strvec[i]==strvec[j]){
                        IsDuplicate=true;
                        break;
                    }
                }
            }
        }
        if(IsDuplicate){
            return false;
        }
        else{
            return true;
        }
    }
    return true;
}
bool StringUtils::IsCommentLine(string &instr){
    string str;
    str=RemoveStrSpace(instr);
    if(str.compare(0,2,"//")==0){
        return true;
    }
    else{
        return false;
    }
}
//**********************************************************
bool StringUtils::IsBracketMatch(ifstream &in,const int &linenum0){
    // the instr="[xxx]"
    vector<int> header,end;
    header.clear();end.clear();
    header.push_back(1);// instr already contains the first header line
    bool IsMatch=false;
    string str;
    int linenum=linenum0;
    IsMatch=false;
    while(!in.eof()){
        getline(in,str);linenum+=1;
        if(!IsCommentLine(str)){
            if(str.find("[end]")!=string::npos||
               str.find("[END]")!=string::npos){
                end.push_back(2);
            }
            else if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
                    (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
                header.push_back(1);
            }
        }

        if((header.size()==end.size())&&(str.length()<1||IsCommentLine(str)||in.eof())){
            IsMatch=true;
            break;
        }
    }
    GoToLine(in,linenum0);
    return IsMatch;
}
bool StringUtils::IsBracketMatch(ifstream &in,const int &linenum0,int &lastend_linenum){
    // the instr="[xxx]"
    vector<int> header,end;
    header.clear();end.clear();
    header.push_back(1);// instr already contains the first header line
    bool IsMatch=false;
    string str;
    int linenum=linenum0;
    IsMatch=false;
    while(!in.eof()){
        getline(in,str);linenum+=1;
        if(!IsCommentLine(str)){
            if(str.find("[end]")!=string::npos||
               str.find("[END]")!=string::npos){
                end.push_back(2);
                lastend_linenum=linenum;
            }
            else if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
                    (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
                header.push_back(1);
            }
        }

        if((header.size()==end.size())&&(str.length()<1||IsCommentLine(str)||in.eof())){
            IsMatch=true;
            break;
        }
    }
    GoToLine(in,linenum0);
    return IsMatch;
}
//****************************************************************
vector<double> StringUtils::SplitStrNum(string &instr){
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
vector<double> StringUtils::SplitStrNum(string &instr,char symbol){
    // cout<<"instr="<<instr<<endl;
    string tempstr;
    tempstr=instr;
    if(symbol!=' '){
        tempstr=RemoveStrSpace(instr);
    }
    else{
        tempstr=instr;
    }
    vector<string> strlist=SplitStr(tempstr,symbol);
    if(strlist.size()>1){
        if(strlist[strlist.size()-1].size()==1){
            if(strlist[strlist.size()-1].at(0)<'0'||strlist[strlist.size()-1].at(0)>'9'){
                strlist.pop_back();
            }
        }
    }
    vector<double> numbers;
    numbers.clear();
    for(int i=0;i<static_cast<int>(strlist.size());i++){
        numbers.push_back(stod(strlist[i]));
    }
    return numbers;
}
vector<double> StringUtils::SplitStrNumAfter(string instr,int pos){
    string str;
    str=instr.substr(pos,instr.size()-pos+1);
    return SplitStrNum(str);
}
//********************************************************
void StringUtils::GoToLine(ifstream &in,const int &linenum){
    in.clear();
    in.seekg(0,ios::beg);// go to the header
    string str;
    for(int i=0;i<linenum;i++) getline(in,str);
}
// for time dependent dirichlet bc
bool StringUtils::IsValidExpression(string str){
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
    if(SplitStrNum(str1).size()>0&&str2.compare(0,2,"*t")==0){
        return true;
    }
    // for t*number case
    str1=str.substr(0,i+1);
    str2=str.substr(i+1);
    if(str1.compare(0,2,"t*")==0&&SplitStrNum(str2).size()>0){
        return true;
    }
    return false;
}