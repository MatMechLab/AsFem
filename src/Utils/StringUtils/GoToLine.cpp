//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/StringUtils.h"

void GoToLine(ifstream &in,const int &linenum)
{
    in.clear();
    in.seekg(0,ios::beg);// go to the header
    string str;
    for(int i=0;i<linenum;i++) getline(in,str);
}