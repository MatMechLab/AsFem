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
//+++ Date   : 2020.06.28
//+++ Purpose: Implement a general message printer for AsFem
//+++          This class is general to provide almost all the
//+++          message print in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/MessagePrinter.h"

MessagePrinter::MessagePrinter(){
}

void MessagePrinter::AsFem_Exit(){
    PrintStars(MessageColor::RED);
    PrintTxt("AsFem exit due to some errors",MessageColor::RED);
    PrintStars(MessageColor::RED);
    PetscEnd();
}

void MessagePrinter::SetColor(const MessageColor &color){
    switch (color) {
        case MessageColor::WHITE:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;37m");// set color to white
            break;
        case MessageColor::RED:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;91m");// set color to bright red
            break;
        case MessageColor::BLUE:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;94m");// set color to bright blue
            break;
        case MessageColor::GREEN:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;32m");// set color to green
            break;
        case MessageColor::YELLOW:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;33m");// set color to white
            break;
        case MessageColor::MAGENTA:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;35m");// set color to white
            break;
        case MessageColor::CYAN:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;36m");// set color to white
            break;
        default:
            break;
    }
}

void MessagePrinter::PrintDashLine(MessageColor color){
    SetColor(color);
    PetscPrintf(PETSC_COMM_WORLD,"***");
    for(int i=0;i<_nWords-6;i++){
        PetscPrintf(PETSC_COMM_WORLD,"-");
    }
    PetscPrintf(PETSC_COMM_WORLD,"***\n");
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//*********************************************
void MessagePrinter::PrintStars(MessageColor color){
    SetColor(color);
    for(int i=0;i<_nWords;i++){
        PetscPrintf(PETSC_COMM_WORLD,"*");
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//*********************************************
void MessagePrinter::PrintTxt(string str,MessageColor color){
    SetColor(color);
    string _Head="*** ";
    string _End=" !!! ***";
    if(str.length()<=_nWords-_Head.length()-_End.length()){
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());

        int i1=static_cast<int>(_Head.size());
        int i2=static_cast<int>(_End.size());
        int i3=static_cast<int>(str.size());

        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    }
    else{
        string substr1,substr2;
        substr1=str.substr(0,_nWords-_Head.length()-_End.length());
        substr2=str.substr(_nWords-_Head.length()-_End.length());
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",substr1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());

        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",substr2.c_str());
        int i1=static_cast<int>(_Head.size());
        int i2=static_cast<int>(_End.size());
        int i3=static_cast<int>(substr2.size());
        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    }
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//**********************************************************
void MessagePrinter::PrintShortTxt(string str,MessageColor color){
    SetColor(color);
    string _Head="*** ";
    string _End=" !!! ***";
    if(str.length()<=_nWords-_Head.length()-_End.length()){
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());

        int i1=static_cast<int>(_Head.size());
        int i2=static_cast<int>(_End.size());
        int i3=static_cast<int>(str.size());

        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    }
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//**********************************************************
void MessagePrinter::PrintLongTxt(string str,MessageColor color){
    SetColor(color);
    MessagePrinter printer;
    vector<string> strvec;
    strvec=printer.SplitStr2Vec(str);
    for(const auto &it:strvec){
        PrintShortTxt(it,color);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//**********************************************************
void MessagePrinter::PrintWelcomeTxt(string str){
    SetColor(MessageColor::CYAN);

    string _Head="*** ";
    string _End =" ***";
    
    PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
    PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());

    int i1=static_cast<int>(_Head.size());
    int i2=static_cast<int>(_End.size());
    int i3=static_cast<int>(str.size());

    for(int i=0;i<_nWords-i1-i2-i3;i++){
        PetscPrintf(PETSC_COMM_WORLD," ");
    }
    PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());

    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//****************************************
vector<string> MessagePrinter::SplitStr2Vec(string str){
    string _Head="*** ";
    string _End=" !!! ***";
    int i1=static_cast<int>(_Head.size());
    int i2=static_cast<int>(_End.size());
    int i3=static_cast<int>(str.size());
    vector<string> strvec;
    strvec.clear();
    if(i3<=_nWords-i1-i2){
        strvec.push_back(str);
    }
    else{
        strvec.clear();
        string substr;
        substr.clear();
        int count=0;
        int nWords=_nWords-i1-i2;
        for(int i=0;i<static_cast<int>(str.length());i++){
            substr.push_back(str.at(i));
            count+=1;
            if(static_cast<int>(substr.size())==nWords){
                strvec.push_back(substr);
                substr.clear();
                if(count==static_cast<int>(str.length())-1){
                    break;
                }
            }
            else{
                if(count==static_cast<int>(str.length())-1){
                    strvec.push_back(substr);
                    break;
                }
            }
        }
    }
    return strvec;
}

//**************************************************
//*** for error message print
//**************************************************
vector<string> MessagePrinter::SplitErrorStr2Vec(string str){
    // string _Head="*** Error:";
    string _Head="***       ";
    string _End=" !!! ***";
    int i1=static_cast<int>(_Head.size());
    int i2=static_cast<int>(_End.size());
    int i3=static_cast<int>(str.size());
    vector<string> strvec;
    strvec.clear();
    if(i3<=_nWords-i1-i2){
        strvec.push_back(str);
    }
    else{
        strvec.clear();
        string substr;
        substr.clear();
        int count=0;
        int nWords=_nWords-i1-i2;
        for(int i=0;i<static_cast<int>(str.length());i++){
            substr.push_back(str.at(i));
            count+=1;
            if(static_cast<int>(substr.size())==nWords){
                strvec.push_back(substr);
                substr.clear();
            }
            else{
                if(count==static_cast<int>(str.length())){
                    strvec.push_back(substr);
                    substr.clear();
                    break;
                }
            }
            if(count==static_cast<int>(str.length())){
                break;
            }
        }
    }
    return strvec;
}
//****************************************************
void MessagePrinter::PrintErrorTxt(string str,bool flag){
    SetColor(MessageColor::RED);
    string _Head ="***       ";
    string _Head1="*** Error:";
    string _End=" !!! ***";
    int i1=static_cast<int>(_Head1.size());
    int i2=static_cast<int>(_End.size());
    int i3=static_cast<int>(str.size());
    if(i3<=_nWords-i1-i2){
        if(flag) PrintStars(MessageColor::RED);
        SetColor(MessageColor::RED);
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());
        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
        if(flag) PrintStars(MessageColor::RED);
        SetColor(MessageColor::RED);
    }
    else{
        string substr1,substr2;
        vector<string> strvec;
        substr1=str.substr(0,_nWords-i1-i2);
        substr2=str.substr(_nWords-i1-i2);
        // for the first line
        if(flag) PrintStars(MessageColor::RED);
        SetColor(MessageColor::RED);
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",substr1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
        // for the later lines
        MessagePrinter printer;
        strvec=printer.SplitErrorStr2Vec(substr2);
        i1=static_cast<int>(_Head.size());
        for(const auto &it:strvec){
            // cout<<"str("<<it.size()<<")="<<it<<endl;
            PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
            PetscPrintf(PETSC_COMM_WORLD,"%s",it.c_str());
            i3=static_cast<int>(it.size());
            for(int i=0;i<_nWords-i1-i2-i3;i++){
                PetscPrintf(PETSC_COMM_WORLD," ");
            }
            PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
        }
        if(flag) PrintStars(MessageColor::RED);
        SetColor(MessageColor::RED);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//************************************************
void MessagePrinter::PrintWarningTxt(string str,bool flag){
    SetColor(MessageColor::YELLOW);
    string _Head ="***         ";
    string _Head1="*** Warning:";
    string _End=" !!! ***";
    int i1=static_cast<int>(_Head1.size());
    int i2=static_cast<int>(_End.size());
    int i3=static_cast<int>(str.size());
    if(i3<=_nWords-i1-i2){
        if(flag) PrintStars(MessageColor::YELLOW);
        SetColor(MessageColor::YELLOW);
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());
        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
        if(flag) PrintStars(MessageColor::YELLOW);
        SetColor(MessageColor::YELLOW);
    }
    else{
        string substr1,substr2;
        vector<string> strvec;
        substr1=str.substr(0,_nWords-i1-i2);
        substr2=str.substr(_nWords-i1-i2);
        // for the first line
        if(flag) PrintStars(MessageColor::YELLOW);
        SetColor(MessageColor::YELLOW);
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",substr1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
        // for the later lines
        MessagePrinter printer;
        strvec=printer.SplitErrorStr2Vec(substr2);
        i1=static_cast<int>(_Head.size());
        for(const auto &it:strvec){
            // cout<<"str("<<it.size()<<")="<<it<<endl;
            PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
            PetscPrintf(PETSC_COMM_WORLD,"%s",it.c_str());
            i3=static_cast<int>(it.size());
            for(int i=0;i<_nWords-i1-i2-i3;i++){
                PetscPrintf(PETSC_COMM_WORLD," ");
            }
            PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
        }
        if(flag) PrintStars(MessageColor::YELLOW);
        SetColor(MessageColor::YELLOW);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//**********************************************
vector<string> MessagePrinter::SplitNormalStr2Vec(string str){
    // string _Head="*** Error:";
    string _Head="*** ";
    string _End =" ***";
    int i1=static_cast<int>(_Head.size());
    int i2=static_cast<int>(_End.size());
    int i3=static_cast<int>(str.size());
    vector<string> strvec;
    strvec.clear();
    if(i3<=_nWords-i1-i2){
        strvec.push_back(str);
    }
    else{
        strvec.clear();
        string substr;
        substr.clear();
        int count=0;
        int nWords=_nWords-i1-i2;
        for(int i=0;i<static_cast<int>(str.length());i++){
            substr.push_back(str.at(i));
            count+=1;
            if(static_cast<int>(substr.size())==nWords){
                strvec.push_back(substr);
                substr.clear();
                if(count==static_cast<int>(str.length())-1){
                    break;
                }
            }
            else{
                if(count==static_cast<int>(str.length())-1){
                    strvec.push_back(substr);
                    break;
                }
            }
        }
    }
    return strvec;
}
void MessagePrinter::PrintNormalTxt(string str,MessageColor color){
    SetColor(color);
    string _Head="*** ";
    string _End =" ***";
    int i1=static_cast<int>(_Head.size());
    int i2=static_cast<int>(_End.size());
    int i3=static_cast<int>(str.size());
    if(i3<=_nWords-i1-i2){
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());
        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    }
    else{
        string substr1,substr2;
        vector<string> strvec;
        substr1=str.substr(0,_nWords-i1-i2);
        substr2=str.substr(_nWords-i1-i2);
        // for the first line
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",substr1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
        // for the later lines
        MessagePrinter printer;
        strvec=printer.SplitNormalStr2Vec(substr2);
        i1=static_cast<int>(_Head.size());
        for(const auto &it:strvec){
            // cout<<"str("<<it.size()<<")="<<it<<endl;
            PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
            PetscPrintf(PETSC_COMM_WORLD,"%s",it.c_str());
            i3=static_cast<int>(it.size());
            for(int i=0;i<_nWords-i1-i2-i3;i++){
                PetscPrintf(PETSC_COMM_WORLD," ");
            }
            PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
        }
    }
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//******************************************
void MessagePrinter::PrintErrorInLineNumber(const int &linenumber){
    SetColor(MessageColor::RED);
    string _Head="*** Error:";
    string _End =" ***";
    char buff[35];
    string str;
    int i1=static_cast<int>(_Head.size());
    int i2=static_cast<int>(_End.size());
    snprintf(buff,35," error detected in line-%d",linenumber);
    str=buff;
    int i3=static_cast<int>(str.size());
    PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
    PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());
    for(int i=0;i<_nWords-i1-i2-i3;i++){
        PetscPrintf(PETSC_COMM_WORLD," ");
    }
    PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}