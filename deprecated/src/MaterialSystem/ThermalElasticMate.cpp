#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::ThermalElasticMate(const int &nDim,const double &t,const double &dt,
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

    if(InputParams.size()<11){
        cout<<"*** Error: for thermech mate,eleven params are required!!! ***"<<endl;
        cout<<"***        you should give D,Cref,Omega,E,nu, eigen strain!***"<<endl;
        cout<<"***        eigen strain is 6 components in voigt notation!!***"<<endl;
        Msg_AsFem_Exit();
    }
    // InputParams[0]--->E
    // InputParams[1]--->nu
    double D=InputParams[0];
    double Cref=InputParams[1];
    double Omega=InputParams[2];
    double E=InputParams[3];
    double nu=InputParams[4];
    double v1=InputParams[5];
    double v2=InputParams[6];
    double v3=InputParams[7];
    double v23=InputParams[8];
    double v13=InputParams[9];
    double v12=InputParams[10];

    // for strain
    _Rank2TensorMateValues[0].SetToZeros();
    // for stress
    _Rank2TensorMateValues[1].SetToZeros();
    // for dSigma/dC
    _Rank2TensorMateValues[2].SetToZeros();
    // for elasticity tensor
    _Rank4TensorMateValues[0].SetToZeros();
    _Rank4TensorMateValues[0].SetFromEandNu(E,nu);

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
    
    // for Identitiy rank-2 tensor
    RankTwoTensor I(0.0);
    RankTwoTensor Ic(0.0);
    Ic.SetFromVoigt(v1,v2,v3,v23,v13,v12);
    I.SetToIdentity();
    

    
    RankTwoTensor totalstrain=(_Rank2TensorMateValues[2]+_Rank2TensorMateValues[2].Transpose())*0.5;
    RankTwoTensor thermstrain(0.0);
    if(nDim==2){
        thermstrain=(1.0/3.0)*Omega*(gpU[3-1]-Cref)*Ic;
    }
    else if(nDim==3){
        thermstrain=(1.0/3.0)*Omega*(gpU[4-1]-Cref)*Ic;
    }
    // this is our strain
    _Rank2TensorMateValues[0]=totalstrain-thermstrain;
    
    // this is our stress, which is extremely simple: sigma=C*strain, done!
    _Rank2TensorMateValues[1]=_Rank4TensorMateValues[0]*_Rank2TensorMateValues[0];
    // this is the dStress/dC rank-2 tensor
    _Rank2TensorMateValues[2]=(-1.0/3.0)*Omega*(_Rank4TensorMateValues[0]*I);

    // use the fourth one to calculate the vonMises stress
    _Rank2TensorMateValues[3]=_Rank2TensorMateValues[1];
    // stress deviator tensor
    _Rank2TensorMateValues[3](1,1)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[3](2,2)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[3](3,3)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();


    _Rank2TensorMateValues[4].SetToIdentity();

    _MateValues[0]=D;
    _MateValues[1]=Omega;
    _MateValues[2]=-(1.0/9.0)*Omega*(_Rank4TensorMateValues[0]*I).Trace();
    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    // sij is the deviator tensor of stress
    _MateValues[3]=sqrt(1.5*_Rank2TensorMateValues[3].DoubleDot(_Rank2TensorMateValues[3]));
    // hydrostatic stress
    _MateValues[4]=_Rank2TensorMateValues[1].Trace()/3.0;
    if(nDim==2){
        _MateValues[5]=_Rank2TensorMateValues[1](1,1);//sigxx
        _MateValues[6]=_Rank2TensorMateValues[1](2,2);//sigyy
        _MateValues[7]=_Rank2TensorMateValues[1](1,2);//sigxy

        _MateValues[8]=_Rank2TensorMateValues[0](1,1);//strxx
        _MateValues[9]=_Rank2TensorMateValues[0](2,2);//stryy
        _MateValues[10]=_Rank2TensorMateValues[0](1,2);//strxy
    }
    else if(nDim==3){
        _MateValues[5]=_Rank2TensorMateValues[1](1,1);//sigxx
        _MateValues[6]=_Rank2TensorMateValues[1](2,2);//sigyy
        _MateValues[7]=_Rank2TensorMateValues[1](3,3);//sigzz
        _MateValues[8]=_Rank2TensorMateValues[1](2,3);//sigyz
        _MateValues[9]=_Rank2TensorMateValues[1](1,3);//sigzx
        _MateValues[10]=_Rank2TensorMateValues[1](1,2);//sigxy

        _MateValues[11]=_Rank2TensorMateValues[0](1,1);//strxx
        _MateValues[12]=_Rank2TensorMateValues[0](2,2);//stryy
        _MateValues[13]=_Rank2TensorMateValues[0](3,3);//strzz
        _MateValues[14]=_Rank2TensorMateValues[0](2,3);//stryz
        _MateValues[15]=_Rank2TensorMateValues[0](1,3);//strzx
        _MateValues[16]=_Rank2TensorMateValues[0](1,2);//strxy
    }
}