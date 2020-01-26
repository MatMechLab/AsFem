#include "ICs/ICSystem.h"

void ICSystem::ApplyConstIC(const vector<double> &Params,const int &DofIndex,
    Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U){
    if(Params.size()<1){
        cout<<"*** Error: for ConstIC, you need at least one param!!!***"<<endl;
        cout<<"***        values=constvalue is expected in [ics]!!!  ***"<<endl;
        Msg_AsFem_Exit();
    }
    int iInd;
    for(int i=1;i<=mesh.GetNodesNum();++i){
        iInd=dofHandler.GetIthNodeJthDofIndex(i,DofIndex)-1;
        U.coeffRef(iInd)=Params[0];
    }
}