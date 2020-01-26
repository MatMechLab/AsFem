#include "ICs/ICSystem.h"

void ICSystem::PrintICBlockInfo() const
{
    for(unsigned int i=0;i<_ICBlockList.size();++i)
    {
        _ICBlockList[i].PrintICBlockInfo();
    }
    if(_ICBlockList.size()>0)
    {
        cout<<"***--------------------------------------------------------***"<<endl;
    }
}