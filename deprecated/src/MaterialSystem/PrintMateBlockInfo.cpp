#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::PrintMateBlockInfo() const
{
    for(unsigned int i=0;i<_MaterialBlockList.size();++i)
    {
        _MaterialBlockList[i].PrintMateBlockInfo();
    }
    if(_MaterialBlockList.size()>0)
    {
        cout<<"***--------------------------------------------------------***"<<endl;
    }
}