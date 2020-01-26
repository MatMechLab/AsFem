#include "BCs/BCSystem.h"
#include "MessagePrint/MessagePrint.h"
void BCSystem::ApplyBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const double &t,const double (&ctan)[2],
                Eigen::SparseMatrix<double,Eigen::RowMajor> &K,
                Eigen::VectorXd &RHS,
                Eigen::VectorXd &U){
    for(auto &it:_BCBlockList){
        _DofIndex=dofHandler.GetDofIndexViaName(it._DofName);
        _bcname=it._BoundaryName;
        _bcvalue=it._value;
        if(it._BCElmtName=="preset"){
            ApplyPresetBC(mesh,dofHandler,fe,ctan,_DofIndex,_bcvalue,_bcname,K,RHS,U);
        }
        else if(it._BCElmtName=="neumann"){
            if(it._IsTimeDependent){
                ApplyNeumannBC(mesh,dofHandler,fe,_DofIndex,_bcvalue*t,_bcname,RHS);
            }
            else{
                ApplyNeumannBC(mesh,dofHandler,fe,_DofIndex,_bcvalue,_bcname,RHS);
            }
        }
        else if(it._BCElmtName=="dirichlet"){
            if(it._IsTimeDependent){
                // cout<<"bcvalue="<<_bcvalue<<", time="<<t<<endl;
                ApplyDirichletBC(mesh,dofHandler,_DofIndex,_bcvalue*t,_bcname,K,RHS,U);
            }
            else{
                ApplyDirichletBC(mesh,dofHandler,_DofIndex,_bcvalue,_bcname,K,RHS,U);
            }
        }
        else{
            cout<<"*** Error: unsupported boundary type !!!                   ***"<<endl;
            Msg_AsFem_Exit();
        }
    }
}