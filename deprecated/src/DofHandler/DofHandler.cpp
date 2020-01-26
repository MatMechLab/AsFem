#include "DofHandler/DofHandler.h"

DofHandler::DofHandler()
{
    _nDofs=0;_nNodes=0;_nElmts=0;_nDofsPerNode=0;_nDim=0;
    _DofNameList.clear();
    _DofIndexList.clear();
    _ElmtalDofActiveFlags.clear();
}

//******************************
void DofHandler::PrintDofsInfo() const
{
    cout<<"*** Dofs name=";
    for(unsigned int i=0;i<_DofNameList.size();i++) cout<<_DofNameList[i]<<" ";
    cout<<endl;
    printf("*** Number of dofs = %8d , each node has %2d dofs ***\n",GetDofsNum(),GetDofsNumPerNode());

}