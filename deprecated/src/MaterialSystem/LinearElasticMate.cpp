#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::LinearElasticMate(const int &nDim,const double &t,const double &dt,
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
        cout<<"*** Error: for mechanics mate,two params are required!!!   ***"<<endl;
        cout<<"***        you should give E and nu !!!                    ***"<<endl;
        Msg_AsFem_Exit();
    }
    // InputParams[0]--->E
    // InputParams[1]--->nu

    // for strain
    _Rank2TensorMateValues[0].SetToZeros();
    // for stress
    _Rank2TensorMateValues[1].SetToZeros();
    // for elasticity tensor
    _Rank4TensorMateValues[0].SetToZeros();
    _Rank4TensorMateValues[0].SetFromEandNu(InputParams[0],InputParams[1]);

    // use the third rank2 tensor as template variale to avoid new memory allocate
    // tips: for nonlinear case, please see our demos,
    // where the powerful tensor calculation is used to make the code extremely simple!!!
    _Rank2TensorMateValues[2].SetToZeros();
    if(nDim==2){
        _Rank2TensorMateValues[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2TensorMateValues[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    
    // this is our strain
    _Rank2TensorMateValues[0]=(_Rank2TensorMateValues[2]+_Rank2TensorMateValues[2].Transpose())*0.5;
    
    // this is our stress, which is extremely simple: sigma=C*strain, done!
    _Rank2TensorMateValues[1]=_Rank4TensorMateValues[0]*_Rank2TensorMateValues[0];

    // use the fourth one to calculate the vonMises stress
    _Rank2TensorMateValues[3]=_Rank2TensorMateValues[1];
    // stress deviator tensor
    _Rank2TensorMateValues[3](1,1)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[3](2,2)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[3](3,3)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();


    _Rank2TensorMateValues[4].SetToIdentity();//used in mechanics elements

    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    // sij is the deviator tensor of stress
    _MateValues[0]=sqrt(1.5*_Rank2TensorMateValues[3].DoubleDot(_Rank2TensorMateValues[3]));
    // hydrostatic stress
    _MateValues[1]=_Rank2TensorMateValues[1].Trace()/3.0;
    if(nDim==2){
        _MateValues[2]=_Rank2TensorMateValues[1](1,1);//sigxx
        _MateValues[3]=_Rank2TensorMateValues[1](2,2);//sigyy
        _MateValues[4]=_Rank2TensorMateValues[1](1,2);//sigxy

        _MateValues[5]=_Rank2TensorMateValues[0](1,1);//strxx
        _MateValues[6]=_Rank2TensorMateValues[0](2,2);//stryy
        _MateValues[7]=_Rank2TensorMateValues[0](1,2);//strxy
    }
    else if(nDim==3){
        _MateValues[2]=_Rank2TensorMateValues[1](1,1);//sigxx
        _MateValues[3]=_Rank2TensorMateValues[1](2,2);//sigyy
        _MateValues[4]=_Rank2TensorMateValues[1](3,3);//sigzz
        _MateValues[5]=_Rank2TensorMateValues[1](2,3);//sigyz
        _MateValues[6]=_Rank2TensorMateValues[1](1,3);//sigzx
        _MateValues[7]=_Rank2TensorMateValues[1](1,2);//sigxy

        _MateValues[ 8]=_Rank2TensorMateValues[0](1,1);//strxx
        _MateValues[ 9]=_Rank2TensorMateValues[0](2,2);//stryy
        _MateValues[10]=_Rank2TensorMateValues[0](3,3);//strzz
        _MateValues[11]=_Rank2TensorMateValues[0](2,3);//stryz
        _MateValues[12]=_Rank2TensorMateValues[0](1,3);//strzx
        _MateValues[13]=_Rank2TensorMateValues[0](1,2);//strxy
    }
}