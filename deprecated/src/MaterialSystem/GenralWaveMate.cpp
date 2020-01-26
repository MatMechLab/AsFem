#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::GenralWaveMate(const int &nDim,const double &t,const double &dt,
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
        cout<<"*** Error: for wave mate,two params are required!!!        ***"<<endl;
        cout<<"***        v1-->mate mode;v2-->wave velocity !!!           ***"<<endl;
        Msg_AsFem_Exit();
    }
    // InputParams[0]--->mate mode
    // InputParams[1]--->wave propagation speed
    
    
   

    if(int(InputParams[0])==1){
        // for mode-1, we use constant wave material
        _MateValues[0]=InputParams[1];// M=D*u*(1-c)
        _MateValues[1]=0.0;// dM/du
        _MateValues[2]=0.0;// F
        _MateValues[3]=0.0;// dF/du
    }
    else if(int(InputParams[0])==2){
        // for mode-2, we use nonlinear wave material
        _MateValues[0]=InputParams[1]*(20.0+sin(gpU[0]));// M=D*u*(1-c)
        _MateValues[1]=InputParams[1]*cos(gpU[0]);// dM/du
        _MateValues[2]=1.0;// F
        _MateValues[3]=0.0;// dF/du
    }
    else if(int(InputParams[0])==3){
        if(InputParams.size()<6){
            cout<<"*** Error: for wave mate-3,six params are required!!!      ***"<<endl;
            cout<<"***        v1-->mate mode;v2-->wave velocity !!!           ***"<<endl;
            cout<<"***        v3-->x0;v4-->y0                   !!!           ***"<<endl;
            cout<<"***        v5-->radius;v6-->applied time     !!!           ***"<<endl;
            Msg_AsFem_Exit();
        }
        // for mode-3, we use central point as source term for a specific time,
        // after that, source term vanish
        _MateValues[0]=InputParams[1];// M=D*u*(1-c)
        _MateValues[1]=0.0;// dM/du
        // only consider about the circle domain
        if(t<=InputParams[5]){
            double dist=(gpCoord(0)-InputParams[2])*(gpCoord(0)-InputParams[2])
                       +(gpCoord(1)-InputParams[3])*(gpCoord(1)-InputParams[3]);
            if(sqrt(dist)<=InputParams[4]){
                const double PI=3.14159265358979323846264338327950;
                _MateValues[2]=0.05*sin(4.0*PI*t);// F
            }
            else{
                _MateValues[2]=0.0;// F
            }
        }
        else{
            _MateValues[2]=0.0;// F
        }
         _MateValues[3]=0.0;// dF/du
    }

}