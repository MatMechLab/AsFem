//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/StringUtils.h"

string RemoveStrSpace(string &instr)
{
    string outstr=instr;
    outstr.erase(remove(outstr.begin(),outstr.end(),' '),outstr.end());
    outstr.erase(remove(outstr.begin(),outstr.end(),'\t'),outstr.end());
    return outstr;
}
//*************************************
string RemoveSymbolFromStr(string &instr,char symbol){
    string outstr=instr;
    outstr.erase(remove(outstr.begin(),outstr.end(),symbol),outstr.end());
    return outstr;
}