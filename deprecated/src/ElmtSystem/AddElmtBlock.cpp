#include "ElmtSystem/ElmtSystem.h"
#include "MessagePrint/MessagePrint.h"

bool ElmtSystem::AddElmtBlock(const ElmtBlock &elmtBlock)
{
    if(_ElmtBlockList.size()<1)
    {
        _ElmtBlockList.push_back(elmtBlock);
        _nElmtBlock=int(_ElmtBlockList.size());
        return true;
    }
    else
    {
        bool NotInList=true;
        for(unsigned int i=0;i<_ElmtBlockList.size();i++)
        {
            if(_ElmtBlockList[i]._ElmtBlockName==elmtBlock._ElmtBlockName)
            {
                NotInList=false;
                break;
            }
        }
        if(NotInList)
        {
            ElmtBlock tempblock;

            tempblock._ElmtBlockName=elmtBlock._ElmtBlockName;
            tempblock._ElmtDofNameList=elmtBlock._ElmtDofNameList;
            tempblock._ElmtDomainBlockName=elmtBlock._ElmtDomainBlockName;
            tempblock._ElmtMateName=elmtBlock._ElmtMateName;
            tempblock._ElmtTypeName=elmtBlock._ElmtTypeName;
            _ElmtBlockList.push_back(tempblock);
            _nElmtBlock=int(_ElmtBlockList.size());
            return NotInList;
        }
        else
        {
            printf("*** Error: duplicate [%20s] in [elmts]***\n",elmtBlock._ElmtBlockName.c_str());
            Msg_AsFem_Exit();
        }
    }
    return true;
}