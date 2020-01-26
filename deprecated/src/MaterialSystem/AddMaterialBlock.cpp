#include "MaterialSystem/MaterialSystem.h"
#include "MessagePrint/MessagePrint.h"

bool MaterialSystem::AddMaterialBlock(const MaterialBlock &mateblock)
{
    if(_MaterialBlockList.size()<1)
    {
        _MaterialBlockList.push_back(mateblock);
        _nMaterialBlocks=int(_MaterialBlockList.size());
        
        return true;
    }
    else
    {
        bool NotInList=true;
        for(unsigned int i=0;i<_MaterialBlockList.size();i++)
        {
            if(_MaterialBlockList[i]._MateBlockName==mateblock._MateBlockName)
            {
                NotInList=false;
                break;
            }
        }
        if(NotInList)
        {
            MaterialBlock temp;
            temp._MateBlockName=mateblock._MateBlockName;
            temp._MateTypeName=mateblock._MateTypeName;
            temp._MateParams=mateblock._MateParams;
            
            _MaterialBlockList.push_back(temp);
            _nMaterialBlocks=int(_MaterialBlockList.size());
            return NotInList;
        }
        else
        {
            printf("*** Error: duplicate [%20s] in [mate]!***\n",mateblock._MateBlockName.c_str());
            Msg_AsFem_Exit();
        }
    }
    return true;
} 