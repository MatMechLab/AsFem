#include "EquationSystem/EquationSystem.h"

void EquationSystem::CreateSparsityPatterns(DofHandler &dofHandler){
    if(!_IsInit){
        Msg_EquationSystem_NotInitialized();
        Msg_AsFem_Exit();
    }
    int e,i,j;
    vector<int> elDofs;
    _ZeroCoeffList.clear();
    elDofs.resize(dofHandler.GetMaxDofsNumPerElmt(),0);
    int nDofs;
    for(e=1;e<=dofHandler.GetBulkElmtsNum();++e){
        dofHandler.GetIthElmtalDofIndex(e,elDofs);
        nDofs=dofHandler.GetIthElmtDofsNum(e);
        for(i=0;i<nDofs;++i){
            for(j=0;j<nDofs;++j){
                _ZeroCoeffList.push_back(T(elDofs[i]-1,elDofs[j]-1,0.0));
            }
        }
    }
    _AMATRIX.setFromTriplets(_ZeroCoeffList.begin(),_ZeroCoeffList.end());
    _AMATRIX.finalize();
    _NNZNum=_AMATRIX.nonZeros();
}