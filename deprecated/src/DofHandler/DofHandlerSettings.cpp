#include "DofHandler/DofHandler.h"

void DofHandler::SetDofNameListFromVec(vector<string> &vec)
{
    _DofNameList.clear();
    int k=0;
    for(unsigned int i=0;i<vec.size();i++)
    {
        _DofNameList.push_back(vec[i]);
        k+=1;
        _DofIndexList.push_back(k);
    }
    _nDofsPerNode=int(vec.size());
}