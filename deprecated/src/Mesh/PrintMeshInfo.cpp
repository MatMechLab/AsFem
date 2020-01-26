#include "Mesh/Mesh.h"

void Mesh::PrintMeshInfo() const
{
    cout<<"***--------------------------------------------------------***"<<endl;
    cout<<"*** Mesh information summary:                              ***"<<endl;
    printf("***   nNodes=%8d,elmt=%8d,bulkelmt=%8d      ***\n",GetNodesNum(),
                                                               GetElmtsNum(),
                                                               GetBulkElmtsNum());
    printf("***   nPhysics=%3d, dim=%3d, min dim=%3d, type=%6s      ***\n",GetPhyGroupsNum(),
                                                                      GetDim(),
                                                                      GetDimMin(),
                                                                      GetMeshType().c_str());
    for(auto it:_PhyIDToNameList)
    {
        printf("***   ID=%6d,name=%16s,nElmts=%8d      ***\n",it.first,it.second.c_str(),
                                                            GetElmtsNumViaPhyName(it.second));
    }
    cout<<"***--------------------------------------------------------***"<<endl;
}
//*******************************************
void Mesh::PrintMeshDetailInfo() const{
    cout<<"***--------------------------------------------------------***"<<endl;
    cout<<"*** Mesh information summary:                              ***"<<endl;
    printf("***   nNodes=%8d,elmt=%8d,bulkelmt=%8d      ***\n",GetNodesNum(),
                                                               GetElmtsNum(),
                                                               GetBulkElmtsNum());
    printf("***   nPhysics=%3d, dim=%3d, min dim=%3d, type=%6s      ***\n",GetPhyGroupsNum(),
                                                                      GetDim(),
                                                                      GetDimMin(),
                                                                      GetMeshType().c_str());
    for(auto it:_PhyIDToNameList)
    {
        printf("***   ID=%6d,name=%16s,nElmts=%8d      ***\n",it.first,it.second.c_str(),
                                                            GetElmtsNumViaPhyName(it.second));
    }
    cout<<"***--------------------------------------------------------***"<<endl;
    for(auto it:_PhyNameToElmtIndexList){
        printf("***   physical name=%16s                       ***\n",it.first.c_str());
        for(int e=1;e<=int(it.second.size());++e){
            int ei=it.second[e-1];
            printf("***     e=%6d:",ei);
            for(int i=1;i<=GetIthElmtNodesNum(ei);++i){
                
                printf("%6d ",GetIthElmtJthConn(ei,i));
                if(i%8==0){
                    printf("\n***              ");
                }
            }
            printf("\n");
        }
    }
    cout<<"***   Nodal coordinate:                                    ***"<<endl;
    for(int i=1;i<=GetNodesNum();++i){
        printf("***    i=%6d: %12.5e  %12.5e  %12.5e  ***\n",i,GetIthNodeJthCoord(i,1),GetIthNodeJthCoord(i,2),GetIthNodeJthCoord(i,3));
    }
    cout<<"***--------------------------------------------------------***"<<endl;
}