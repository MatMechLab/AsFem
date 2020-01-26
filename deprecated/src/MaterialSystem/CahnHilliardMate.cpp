#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::CahnHilliardMate(const int &nDim,const double &t,const double &dt,
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

    if(InputParams.size()<3){
        cout<<"*** Error: for cahnhilliard mate,three params are required!***"<<endl;
        Msg_AsFem_Exit();
    }
    // InputParams[0]--->chi
    // InputParams[1]--->kappa
    
    
   

    double c=gpU[0];// if you use cahn-hilliard,
                    // must make sure the first dof is concentration, not mu!!!
    double tol=1.0e-5;
    if(c<tol) c=tol;
    if(c>1.0-tol) c=1.0-tol;
    _MateValues[0]=InputParams[0]*c*(1-c);// M=D*c*(1-c)
    _MateValues[1]=InputParams[0]*(1-2*c);// dM/dc
    _MateValues[2]=InputParams[1];// Chi
    _MateValues[3]=InputParams[2];// Kappa
    _MateValues[4]=c*log(c)+(1-c)*log(1-c)+InputParams[1]*c*(1-c);// free energy
    _MateValues[5]=log(c)-log(1-c)+InputParams[1]*(1-2*c);// mu=df/dc->chemical potential
    _MateValues[6]=1.0/c+1.0/(1.0-c)-InputParams[1]*2.0;//dmu/dc

}