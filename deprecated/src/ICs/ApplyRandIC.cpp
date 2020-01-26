#include "ICs/ICSystem.h"

void ICSystem::ApplyRandIC(const vector<double> &Params,const int &DofIndex,
    Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U){
    if(Params.size()<2){
        cout<<"*** Error: for RandIC, you need at least two params!!!***"<<endl;
        cout<<"***        values=min max is expected in [ics]!!!     ***"<<endl;
        Msg_AsFem_Exit();
    }
    const double minval=Params[0];
    const double maxval=Params[1];
    int iInd;
    srand(time(0));
    for(int i=1;i<=mesh.GetNodesNum();++i){
        iInd=dofHandler.GetIthNodeJthDofIndex(i,DofIndex)-1;
        U.coeffRef(iInd)=minval+(maxval-minval)*rand()/RAND_MAX;
    }
}