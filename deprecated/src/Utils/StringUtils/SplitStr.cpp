#include "Utils/StringUtils.h"

vector<string> SplitStr(string &instr,char symbol)
{
    vector<string> outstr;
    outstr.clear();

    stringstream ss(instr);
    string temp;
    while (getline(ss,temp,symbol))
    {
        outstr.push_back(temp);
    }
    vector<string> outstr0;
    outstr0.clear();
    if(symbol==' ')
    {
        for(unsigned int i=0;i<outstr.size();i++)
        {
            temp=outstr[i];
            temp.erase(remove(temp.begin(),temp.end(),' '),temp.end());
            if(temp.length()>0)
            {
                outstr0.push_back(temp);
            }
        }
    }
    outstr=outstr0;
    return outstr;
}