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