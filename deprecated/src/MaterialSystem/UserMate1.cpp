#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::UserMate1(const int &nDim,const double &t,const double &dt,
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

    if(InputParams.size()<2){
        cout<<"*** Error: for umat1, two parameters are required!!!       ***"<<endl;
        cout<<"***        you should give E and nu !!!                    ***"<<endl;
        Msg_AsFem_Exit();
    }
    
    //************************
    //*** this is an example for user-defined material(used by uel-->nonlinear geometry problem)
    //**** MateVals[0]-->store E
    //**** MateVals[1]-->store nu
    _MateValues[0]=InputParams[0];// Young's modulus
    _MateValues[1]=InputParams[1];// Poisson ratio
    
}