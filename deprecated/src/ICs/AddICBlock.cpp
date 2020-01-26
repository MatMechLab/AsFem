#include "ICs/ICSystem.h"
#include "MessagePrint/MessagePrint.h"

bool ICSystem::AddICBlock(const ICBlock &icblock)
{
    if(_ICBlockList.size()<1)
    {
        _ICBlockList.push_back(icblock);
        _nICBlocks=int(_ICBlockList.size());
        return true;
    }
    else
    {
        bool NotInList=true;
        for(unsigned int i=0;i<_ICBlockList.size();i++)
        {
            if(_ICBlockList[i]._ICBlockName==icblock._ICBlockName)
            {
                NotInList=false;
                break;
            }
        }
        if(NotInList)
        {
            ICBlock temp;
            temp._ICBlockName=icblock._ICBlockName;
            temp._ICElmtName=icblock._ICElmtName;
            temp._BlockName=icblock._BlockName;
            temp._DofName=icblock._DofName;
            temp._value=icblock._value;
            _ICBlockList.push_back(temp);
            _nICBlocks=int(_ICBlockList.size());
            return NotInList;
        }
        else
        {
            printf("*** Error: duplicate [%20s] in [ics]!!***\n",icblock._ICBlockName.c_str());
            Msg_AsFem_Exit();
        }
    }
    return true;
}