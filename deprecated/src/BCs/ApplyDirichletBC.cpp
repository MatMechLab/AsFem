#include "BCs/BCSystem.h"

void BCSystem::ApplyDirichletBC(Mesh &mesh,DofHandler &dofHandler,
                        const int &DofIndex,const double &bcvalue,string bcname,
                        Eigen::SparseMatrix<double,Eigen::RowMajor> &K,
                        Eigen::VectorXd &RHS,
                        Eigen::VectorXd &U){
    int e,ee,i,j,iInd;
    
    
    for(e=1;e<=mesh.GetElmtsNumViaPhyName(bcname);++e){
        ee=mesh.GetElmtIndexViaPhyName(bcname,e);
        // cout<<"bcname="<<bcname<<", bcvalue="<<bcvalue<<":";
        for(i=1;i<=_nNodesPerBCElmt;++i){
            j=mesh.GetIthElmtJthConn(ee,i);
            iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
            // cout<<iInd+1<<" ";
            U.coeffRef(iInd)=bcvalue;
            RHS.coeffRef(iInd)=0.0;
            K.coeffRef(iInd,iInd)=_MaxKMatrixValue;
        }
        // cout<<endl;
    }
}