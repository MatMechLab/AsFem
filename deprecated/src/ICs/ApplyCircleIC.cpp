#include "ICs/ICSystem.h"

void ICSystem::ApplyCircleIC(const vector<double> &Params,const int &DofIndex,
    Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U){
    if(Params.size()<5){
        cout<<"*** Error: for CircleIC, you need at least five params***"<<endl;
        cout<<"***        values=cx cy r v1 v2 are expected in [ics]!***"<<endl;
        Msg_AsFem_Exit();
    }
    int iInd;
    const double cx=Params[0];
    const double cy=Params[1];
    const double r=Params[2];
    const double v1=Params[3];
    const double v2=Params[4];
    double x,y,dist;
    for(int i=1;i<=mesh.GetNodesNum();++i){
        iInd=dofHandler.GetIthNodeJthDofIndex(i,DofIndex)-1;
        x=mesh.GetIthNodeJthCoord(i,1);
        y=mesh.GetIthNodeJthCoord(i,2);
        dist=sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy));
        if(dist<=r){
            U.coeffRef(iInd)=v1;
        }
        else{
            U.coeffRef(iInd)=v2;
        }
    }
}