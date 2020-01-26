#include "Utils/StringUtils.h"

string RemoveStrSpace(string &instr)
{
    string outstr=instr;
    outstr.erase(remove(outstr.begin(),outstr.end(),' '),outstr.end());
    outstr.erase(remove(outstr.begin(),outstr.end(),'\t'),outstr.end());
    return outstr;
}