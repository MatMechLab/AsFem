#include "BCs/BCSystem.h"

void BCSystem::ApplyPreSetDirichletBC(const double &t,Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U){
    int e,ee,i,j,iInd;
    for(auto &it:_BCBlockList){
        _DofIndex=dofHandler.GetDofIndexViaName(it._DofName);
        _bcname=it._BoundaryName;
        _bcvalue=it._value;
        if(it._IsTimeDependent){
            _bcvalue=it._value*t;
        }
        if(it._BCElmtName=="preset"||it._BCElmtName=="dirichlet"){
            for(e=1;e<=mesh.GetElmtsNumViaPhyName(_bcname);++e){
                ee=mesh.GetElmtIndexViaPhyName(_bcname,e);
                // cout<<"bcname="<<_bcname<<", bcvalue="<<_bcvalue<<", e="<<ee<<":";
                for(i=1;i<=_nNodesPerBCElmt;++i){
                    j=mesh.GetIthElmtJthConn(ee,i);
                    iInd=dofHandler.GetIthNodeJthDofIndex(j,_DofIndex)-1;
                    U.coeffRef(iInd)=_bcvalue;
                    // cout<<iInd+1<<" ";
                }
                // cout<<endl;
            }
        }
        else{
            continue;
        }
    }
}