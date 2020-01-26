#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::SaintVenantMate(const int &nDim,const double &t,const double &dt,
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
    const double EE=InputParams[0];
    const double nu=InputParams[1];
    const double lambda=EE*nu/((1+nu)*(1-2*nu));// lame const
    const double mu=EE/(2*(1+nu));

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
    RankTwoTensor F(0.0),I(0.0),C(0.0),Cinv(0.0),E(0.0);
    RankTwoTensor pk2(0.0);
    I.SetToIdentity();
    
    if(nDim==2){
        _Rank2TensorMateValues[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2TensorMateValues[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    F=I+_Rank2TensorMateValues[2];// F=I+U_{i,j}
    C=F.Transpose()*F;//C=F^tF
    Cinv=C.Inverse();
    E=(C-I)*0.5;

    // now our strain is Green-Lagrange finite strain
    _Rank2TensorMateValues[0]=E;
    
    // Second Piola Kirchhoff Stress tensor
    pk2=lambda*E.Trace()*I+2*mu*E;
    // pk2=(0.5*lambda*(Cinv.Trace()-nDim)-mu)*I+Cinv*mu;
    // pk2.Print();
    // cout<<endl;
    // this pk1 stress is used in mechanics element
    _Rank2TensorMateValues[1]=F*pk2;// pk1 stress(since we use reference configuration)
    // _Rank2TensorMateValues[1]=pk2;
    _Rank2TensorMateValues[4]=F;

    // now we calculate the Cijkl for saintVenant
    // Cijkl=lamda*delta_ij*delta_kl+2*mu*0.5*(delta_ik*delta_jl+delta_il*delta_jk)
    _Rank4TensorMateValues[0]=lambda*I.CrossDot(I)+mu*I.ODot(I);
    // _Rank4TensorMateValues[0].SetFromEandNu(EE,nu);
    _Rank2TensorMateValues[1]=_Rank4TensorMateValues[0]*E;

    // use the fourth one to calculate the vonMises stress
    _Rank2TensorMateValues[3]=(_Rank2TensorMateValues[1]*F.Transpose())*(1.0/F.Det());//cauch stress=FPF^T/J
    // stress deviator tensor
    _Rank2TensorMateValues[3](1,1)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[3](2,2)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[3](3,3)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();

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