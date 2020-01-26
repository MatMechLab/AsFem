#include "DofHandler/DofHandler.h"

void DofHandler::PrintDofInfo() const{
    
    cout<<"*** DofHandler information summary:                        ***"<<endl;
    printf("***   nDofs=%7d,ActiveDofs=%7d,nDofsPerNode=%2d     ***\n",GetDofsNum(),
                                                                 GetActiveDofsNum(),
                                                                 GetDofsNumPerNode());
    cout<<"***--------------------------------------------------------***"<<endl;
}

//*************************************
void DofHandler::PrintDofMap() const{
    cout<<"*** DofHandler information summary:                        ***"<<endl;
    cout<<"*** DofMap print:                                          ***"<<endl;
    for(int e=1;e<=GetBulkElmtsNum();++e){
        printf("*** %6d-th elmt:",e);
        for(int i=1;i<=GetIthElmtDofsNum(e);++i){
            printf("%6d ",GetIthElmtJthDofIndex(e,i));
        }
        cout<<endl;
    }
    cout<<"***--------------------------------------------------------***"<<endl;
    
}