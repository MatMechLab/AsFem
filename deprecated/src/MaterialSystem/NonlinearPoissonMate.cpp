#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::NonlinearPoissonMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld){
    
    // just to get rid of warning to error compile mistakes
    // for normal user, you don't need to care about this(warning will always be warning)
    // for developer, please always remember to enable "warning to error" compile flags!!!
    if(t||dt){}
    if(gpCoord(0)){}
    if(gpU[0]){}
    if(gpV[0]){}
    if(gpGradU[0](0)){}
    if(gpGradV[0](0)){}
    if(gpHist[0]){}
    if(gpHistOld[0]){}
    if(InputParams[0]){}
    
    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    //**** MateVals[0]-->store sigma
    //**** MateVals[1]-->store dsigma/dphi(for constant case, it is zero)
    //**** MateVals[2]-->store F
    //**** MateVals[3]-->store dF/dphi (for constant case, it is zero)
    if(nDim==1){
        _MateValues[0]=2.0+sin(gpCoord(0));// sigma=2.0*sin(x)*cos(y)
    }
    else if(nDim==2){
        _MateValues[0]=3.0+sin(gpU[0])*sin(gpCoord(0)*gpCoord(1));// sigma=2.0*sin(x)*cos(y)
        _MateValues[1]=cos(gpU[0])*sin(gpCoord(0)*gpCoord(1));// dsigma/dphi
    }
    else if(nDim==3){
        _MateValues[0]=2.0+cos(gpU[0])*sin(gpCoord(0)*gpCoord(1)*gpCoord(2));// sigma=2.0*sin(x)*cos(y)
        _MateValues[1]=-sin(gpU[0])*sin(gpCoord(0)*gpCoord(1)*gpCoord(2));
    }
    
    _MateValues[2]=2.0+gpU[0]*sin(gpCoord(0))*cos(gpCoord(1));// F
    _MateValues[3]=sin(gpCoord(0))*cos(gpCoord(1));// dF/dphi

    // in this case, the equation becomes:
    // lap(phi)=1
    // so now you can see that, in AsFem, we use this approach
    // to allow modeling of different models
    // one uel+one umate, which is quite helpful for simulations.
}