#include "Utils/StringUtils.h"

string RemoveSymbolFromStr(string &instr,char symbol)
{
    string outstr=instr;
    outstr.erase(remove(outstr.begin(),outstr.end(),symbol),outstr.end());
    return outstr;
}