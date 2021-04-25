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
//+++ Date   : 2020.06.28
//+++ Purpose: Implement a general message printer for AsFem
//+++          This class is general to provide almost all the
//+++          message print in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

#include "petsc.h"

#include "Utils/MessageColor.h"

using namespace std;


class MessagePrinter{
public:
    MessagePrinter();

    static void PrintTxt(string str,MessageColor color=MessageColor::WHITE);
    static void PrintShortTxt(string str,MessageColor color=MessageColor::WHITE);
    static void PrintLongTxt(string str,MessageColor color=MessageColor::WHITE);
    static void PrintErrorTxt(string str,bool flag=true);
    static void PrintWarningTxt(string str,bool flag=true);

    static void PrintNormalTxt(string str,MessageColor color=MessageColor::WHITE);
    
    static void PrintWelcomeTxt(string str);
    static void PrintStars(MessageColor color=MessageColor::WHITE);
    static void PrintDashLine(MessageColor color=MessageColor::WHITE);

    //***********************************
    //*** for input file read
    //***********************************
    static void PrintErrorInLineNumber(const int &linenumber);

    static void AsFem_Exit();

    static void SetColor(const MessageColor &color);

private:
    static const int _nWords=77;
    vector<string> SplitStr2Vec(string str);
    vector<string> SplitErrorStr2Vec(string str);
    vector<string> SplitNormalStr2Vec(string str);

};