#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::ConstDiffusionMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld){
    
    // just to get rid of warning to error compile mistakes
    // for normal user, you don't need to care about this(warning will always be warning)
    // for developer, please always remember to enable "warning to error" compile flags!!!
    if(nDim){}
    if(t||dt){}
    if(gpCoord(0)){}
    if(gpU[0]){}
    if(gpV[0]){}
    if(gpGradU[0](0)){}
    if(gpGradV[0](0)){}
    if(gpHist[0]){}
    if(gpHistOld[0]){}

    if(InputParams.size()<1){
        cout<<"*** Error: for constdiffusion mate, one param is required!!***"<<endl;
        Msg_AsFem_Exit();
    }
    
    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    //**** MateVals[0]-->store D
    //**** MateVals[1]-->store dD/dc(for constant case, it is zero)
    _MateValues[0]=InputParams[0];// D
    _MateValues[1]=0.0;// dD/dc

    // in this case, the equation becomes:
    // lap(phi)=1
    // so now you can see that, in AsFem, we use this approach
    // to allow modeling of different models
    // one uel+one umate, which is quite helpful for simulations.
}