//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/StringUtils.h"

bool IsBracketMatch(ifstream &in,const int &linenum0)
{
    // the instr="[xxx]"
    vector<int> header,end;
    header.clear();end.clear();
    header.push_back(1);// instr already contains the first header line
    bool IsMatch=false;
    string str;
    int linenum=linenum0;
    IsMatch=false;
    while(!in.eof())
    {
        getline(in,str);linenum+=1;
        if(!IsCommentLine(str))
        {
            if(str.find("[end]")!=string::npos||
               str.find("[END]")!=string::npos)
            {
                end.push_back(2);
            }
            else if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
                    (str.find("[end]")==string::npos&&str.find("[END]")==string::npos))
            {
                header.push_back(1);
            }
        }

        if((header.size()==end.size())&&(str.length()<1||IsCommentLine(str)||in.eof()))
        {
            IsMatch=true;
            break;
        }
    }

    GoToLine(in,linenum0);

    return IsMatch;

}

//*************************************
bool IsBracketMatch(ifstream &in,const int &linenum0,int &lastend_linenum)
{
    // the instr="[xxx]"
    vector<int> header,end;
    header.clear();end.clear();
    header.push_back(1);// instr already contains the first header line
    bool IsMatch=false;
    string str;
    int linenum=linenum0;
    IsMatch=false;
    while(!in.eof())
    {
        getline(in,str);linenum+=1;
        if(!IsCommentLine(str))
        {
            if(str.find("[end]")!=string::npos||
               str.find("[END]")!=string::npos)
            {
                end.push_back(2);
                lastend_linenum=linenum;
            }
            else if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
                    (str.find("[end]")==string::npos&&str.find("[END]")==string::npos))
            {
                header.push_back(1);
            }
        }

        if((header.size()==end.size())&&(str.length()<1||IsCommentLine(str)||in.eof()))
        {
            IsMatch=true;
            break;
        }
    }

    GoToLine(in,linenum0);

    return IsMatch;
}