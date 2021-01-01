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
//+++ Date   : 2020.06.30
//+++ Purpose: Implement a general class for string manipulate
//+++          For example, string cases convert and the split..
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <vector>
#include <sstream>

using namespace std;

class StringUtils{
public:
    StringUtils();

    static string StrToLower(string instr);
    static string StrToUpper(string instr);

    static string RemoveStrSpace(string &instr);
    static string RemoveSymbolFromStr(string &instr,char symbol);
    static vector<string> SplitStr(string &instr,char symbol);
    
    static bool IsUniqueStrVec(vector<string> &strvec);
    static bool IsCommentLine(string &instr);
    
    static bool IsBracketMatch(ifstream &in,const int &linenum0);
    static bool IsBracketMatch(ifstream &in,const int &linenum0,int &lastend_linenum);
    
    // for string to number convert
    static vector<double> SplitStrNum(string &instr);
    static vector<double> SplitStrNum(string &instr,char symbol);
    static vector<double> SplitStrNumAfter(string instr,int pos);
    
    static void GoToLine(ifstream &in,const int &linenum);
    // for time dependent dirichlet bc
    static bool IsValidExpression(string str);

};