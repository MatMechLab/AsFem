#include "Utils/StringUtils.h"

void GoToLine(ifstream &in,const int &linenum)
{
    in.clear();
    in.seekg(0,ios::beg);// go to the header
    string str;
    for(int i=0;i<linenum;i++) getline(in,str);
}