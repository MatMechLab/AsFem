//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/StringUtils.h"

bool IsCommentLine(string &instr)
{
    string str;
    str=RemoveStrSpace(instr);
    if(str.compare(0,2,"//")==0)
    {
        return true;
    }
    else
    {
        return false;
    }
}