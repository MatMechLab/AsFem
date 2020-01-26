#include "Utils/StringUtils.h"

vector<double> SplitStrNumAfter(string instr,int pos)
{
    string str;
    str=instr.substr(pos,instr.size()-pos+1);
    return SplitStrNum(str);
}