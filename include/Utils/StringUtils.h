//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_STRINGUTILS_H
#define ASFEM_STRINGUTILS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <vector>
#include <sstream>

using namespace std;

string StrToLower(string instr);
string StrToUpper(string instr);

string RemoveStrSpace(string &instr);
string RemoveSymbolFromStr(string &instr,char symbol);
vector<string> SplitStr(string &instr,char symbol);

bool IsUniqueStrVec(vector<string> &strvec);

bool IsCommentLine(string &instr);

bool IsBracketMatch(ifstream &in,const int &linenum0);
bool IsBracketMatch(ifstream &in,const int &linenum0,int &lastend_linenum);

// for string to number convert
vector<double> SplitStrNum(string &instr);
vector<double> SplitStrNum(string &instr,char symbol);
vector<double> SplitStrNumAfter(string instr,int pos);


void GoToLine(ifstream &in,const int &linenum);

// for time dependent dirichlet bc
bool IsValidExpression(string str);

#endif // ASFEM_STRINGUTILS_H