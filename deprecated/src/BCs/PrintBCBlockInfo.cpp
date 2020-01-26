#include "BCs/BCSystem.h"

void BCSystem::PrintBCBlockInfo() const
{
    for(unsigned int i=0;i<_BCBlockList.size();i++)
    {
        _BCBlockList[i].PrintBCBlockInfo();
    }
    if(_BCBlockList.size()>0)
    {
        cout<<"***--------------------------------------------------------***"<<endl;
    }
}