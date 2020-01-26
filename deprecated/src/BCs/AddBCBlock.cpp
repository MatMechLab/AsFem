#include "BCs/BCSystem.h"
#include "MessagePrint/MessagePrint.h"

bool BCSystem::AddBCBlock(const BCBlock &bcblock)
{
    if(_BCBlockList.size()<1)
    {
        _BCBlockList.push_back(bcblock);
        _nBCBlocks=int(_BCBlockList.size());
        return true;
    }
    else
    {
        bool NotInList=true;
        for(unsigned int i=0;i<_BCBlockList.size();i++)
        {
            if(_BCBlockList[i]._BCBlockName==bcblock._BCBlockName)
            {
                NotInList=false;
                break;
            }
        }
        if(NotInList)
        {
            BCBlock temp;
            temp._BCBlockName=bcblock._BCBlockName;
            temp._BCElmtName=bcblock._BCElmtName;
            temp._BoundaryName=bcblock._BoundaryName;
            temp._DofName=bcblock._DofName;
            temp._value=bcblock._value;
            temp._IsTimeDependent=bcblock._IsTimeDependent;
            _BCBlockList.push_back(temp);
            _nBCBlocks=int(_BCBlockList.size());
            return NotInList;
        }
        else
        {
            printf("*** Error: duplicate [%20s] in [bcs]!!***\n",bcblock._BCBlockName.c_str());
            Msg_AsFem_Exit();
        }
    }
    return true;
}