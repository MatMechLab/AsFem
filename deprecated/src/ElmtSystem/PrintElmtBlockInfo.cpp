#include "ElmtSystem/ElmtSystem.h"


void ElmtSystem::PrintElmtBlockInfo() const
{
    for(unsigned int i=0;i<_ElmtBlockList.size();++i)
    {
        _ElmtBlockList[i].PrintInfo();
    }
    if(_ElmtBlockList.size()>0)
    {
        cout<<"***--------------------------------------------------------***"<<endl;
    }
}