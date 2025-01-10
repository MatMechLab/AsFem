//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
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

void MessagePrinter::exitAsFem(){
    printStars(MessageColor::RED);
    printTxt("AsFem exit due to some errors",MessageColor::RED);
    printStars(MessageColor::RED);
    PetscEnd();
}

void MessagePrinter::setColor(const MessageColor &color){
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
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;33m");// set color to yellow
            break;
        case MessageColor::MAGENTA:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;35m");// set color to megenta
            break;
        case MessageColor::CYAN:
            PetscPrintf(PETSC_COMM_WORLD,"\033[1;36m");// set color to cyan
            break;
        default:
            break;
    }
}

void MessagePrinter::printDashLine(MessageColor color){
    setColor(color);
    PetscPrintf(PETSC_COMM_WORLD,"***");
    for(int i=0;i<_nWords-6;i++){
        PetscPrintf(PETSC_COMM_WORLD,"-");
    }
    PetscPrintf(PETSC_COMM_WORLD,"***\n");
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//*********************************************
void MessagePrinter::printStars(MessageColor color){
    setColor(color);
    for(int i=0;i<_nWords;i++){
        PetscPrintf(PETSC_COMM_WORLD,"*");
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
void MessagePrinter::printSingleTxt(string str,MessageColor color){
    setColor(color);
    PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());
    setColor(MessageColor::WHITE);
}
void MessagePrinter::printNewLine(){
    PetscPrintf(PETSC_COMM_WORLD,"\n");
}
//*********************************************
void MessagePrinter::printTxt(string str,MessageColor color){
    setColor(color);
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
void MessagePrinter::printErrorTxt(string str,bool flag){
    const string HeadTxt="*** Error: ";
    const string EndTxt=" !!! ***";
    vector<string> subvec;
    MessagePrinter printer;
    subvec=printer.splitStr2Vec(HeadTxt,EndTxt,str);
    int i1=static_cast<int>(HeadTxt.size());
    int i2=static_cast<int>(EndTxt.size());
    int i3=static_cast<int>(str.size());
    setColor(MessageColor::RED);
    if(flag){
        printStars(MessageColor::RED);
    }
    for(const auto &it:subvec){
        setColor(MessageColor::RED);
        PetscPrintf(PETSC_COMM_WORLD,"%s",HeadTxt.c_str());
        setColor(MessageColor::WHITE);
        PetscPrintf(PETSC_COMM_WORLD,"%s",it.c_str());
        i3=static_cast<int>(it.size());
        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        setColor(MessageColor::RED);
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",EndTxt.c_str());
    }
    if(flag){
        printStars(MessageColor::RED);
    }
}
//**********************************************************
void MessagePrinter::printWarningTxt(string str,bool flag){
    const string HeadTxt="*** Warning: ";
    const string EndTxt=" !!! ***";
    vector<string> subvec;
    MessagePrinter printer;
    subvec=printer.splitStr2Vec(HeadTxt,EndTxt,str);
    int i1=static_cast<int>(HeadTxt.size());
    int i2=static_cast<int>(EndTxt.size());
    int i3=static_cast<int>(str.size());
    setColor(MessageColor::YELLOW);
    if(flag){
        printStars(MessageColor::YELLOW);
    }
    for(const auto &it:subvec){
        setColor(MessageColor::YELLOW);
        PetscPrintf(PETSC_COMM_WORLD,"%s",HeadTxt.c_str());
        setColor(MessageColor::WHITE);
        PetscPrintf(PETSC_COMM_WORLD,"%s",it.c_str());
        i3=static_cast<int>(it.size());
        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        setColor(MessageColor::YELLOW);
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",EndTxt.c_str());
    }
    if(flag){
        printStars(MessageColor::YELLOW);
    }
}
//**********************************************************
void MessagePrinter::printShortTxt(string str,MessageColor color){
    setColor(color);
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
void MessagePrinter::printLongTxt(string str,MessageColor color){
    setColor(color);
    MessagePrinter printer;
    vector<string> strvec;
    strvec=printer.splitStr2Vec("***"," !!! ***",str);
    for(const auto &it:strvec){
        printShortTxt(it,color);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//**********************************************************
void MessagePrinter::printWelcomeTxt(string str){
    setColor(MessageColor::BLUE);

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
//**************************************
void MessagePrinter::printNormalTxt(string str,MessageColor color){
    setColor(color);
    string HeadTxt="*** ";
    string EndTxt =" ***";
    int i1=static_cast<int>(HeadTxt.size());
    int i2=static_cast<int>(EndTxt.size());
    int i3=static_cast<int>(str.size());

    vector<string> strvec;
    
    MessagePrinter printer;
    strvec=printer.splitStr2Vec(HeadTxt,EndTxt,str);
    for(const auto &it:strvec){
        setColor(color);
        PetscPrintf(PETSC_COMM_WORLD,"%s",HeadTxt.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",it.c_str());
        i3=static_cast<int>(it.size());
        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",EndTxt.c_str());
    }
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//******************************************
void MessagePrinter::printErrorInLineNumber(const int &linenumber){
    setColor(MessageColor::RED);
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
    setColor(MessageColor::WHITE);
    PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());
    for(int i=0;i<_nWords-i1-i2-i3;i++){
        PetscPrintf(PETSC_COMM_WORLD," ");
    }
    setColor(MessageColor::RED);
    PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    PetscPrintf(PETSC_COMM_WORLD,"\033[0m");// recover color
}
//**********************************************
vector<string> MessagePrinter::splitStr2Vec(const string &headtxt,const string &endtxt,string str){
    int i1=static_cast<int>(headtxt.size());
    int i2=static_cast<int>(endtxt.size());
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
                if(count==static_cast<int>(str.length())){
                    break;
                }
            }
            else{
                if(count==static_cast<int>(str.length())){
                    strvec.push_back(substr);
                    break;
                }
            }
        }
    }
    return strvec;
}