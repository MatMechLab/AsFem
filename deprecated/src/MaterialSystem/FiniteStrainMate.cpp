#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::FiniteStrainMate(const int &nDim,const double &t,const double &dt,
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
    RankTwoTensor F(0.0),I(0.0),C(0.0),E(0.0);
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
    E=(C-I)*0.5;
    pk2=_Rank4TensorMateValues[0]*E;

    _Rank2TensorMateValues[4]=F;
    //_Rank2TensorMateValues[4].SetToIdentity();
 

    // now our strain is Green-Lagrange finite strain
    _Rank2TensorMateValues[0]=E;
   
    
    // now we calculate the Cijkl
    _Rank2TensorMateValues[1]=F*pk2;// bad convergence behavior
    //_Rank2TensorMateValues[1]=_Rank4TensorMateValues[0]*E;
    

    if(InputParams.size()==3){
        if(InputParams[2]==1){
            cout<<"displacement="<<gpU[0]<<", "<<gpU[1]<<endl;
            cout<<"gradU="<<endl;
    cout<<gpGradU[0](0)<<","
        <<gpGradU[0](1)<<","
        <<gpGradU[0](2)<<endl
        <<gpGradU[1](0)<<","
        <<gpGradU[1](1)<<","
        <<gpGradU[1](2)<<endl;
    cout<<"displacement gradient:"<<endl;
    cout<<_Rank2TensorMateValues[2](1,1)<<","
        <<_Rank2TensorMateValues[2](1,2)<<","
        <<_Rank2TensorMateValues[2](1,3)<<endl
        <<_Rank2TensorMateValues[2](2,1)<<","
        <<_Rank2TensorMateValues[2](2,2)<<","
        <<_Rank2TensorMateValues[2](2,3)<<endl
        <<_Rank2TensorMateValues[2](3,1)<<","
        <<_Rank2TensorMateValues[2](3,2)<<","
        <<_Rank2TensorMateValues[2](3,3)<<endl<<endl;

    cout<<"F="<<endl;
    cout<<F(1,1)<<","
        <<F(1,2)<<","
        <<F(1,3)<<endl
        <<F(2,1)<<","
        <<F(2,2)<<","
        <<F(2,3)<<endl
        <<F(3,1)<<","
        <<F(3,2)<<","
        <<F(3,3)<<endl<<endl;

    cout<<"E="<<endl;
    cout<<E(1,1)<<","
        <<E(1,2)<<","
        <<E(1,3)<<endl
        <<E(2,1)<<","
        <<E(2,2)<<","
        <<E(2,3)<<endl
        <<E(3,1)<<","
        <<E(3,2)<<","
        <<E(3,3)<<endl<<endl;
    
    cout<<"stress="<<endl;
    cout<<_Rank2TensorMateValues[1](1,1)<<","
        <<_Rank2TensorMateValues[1](1,2)<<","
        <<_Rank2TensorMateValues[1](1,3)<<endl
        <<_Rank2TensorMateValues[1](2,1)<<","
        <<_Rank2TensorMateValues[1](2,2)<<","
        <<_Rank2TensorMateValues[1](2,3)<<endl
        <<_Rank2TensorMateValues[1](3,1)<<","
        <<_Rank2TensorMateValues[1](3,2)<<","
        <<_Rank2TensorMateValues[1](3,3)<<endl<<endl;

    cout<<"Cijkl="<<endl;
    for(int i=1;i<=3;++i){
        for(int j=1;j<=3;++j){
            for(int k=1;k<=3;++k){
                for(int l=1;l<=3;++l){
                    cout<<_Rank4TensorMateValues[0](i,j,k,l)<<",";
                }
            }
            cout<<endl;
        }
    }
    // cout<<_Rank4TensorMateValues[0](1,1,1,1)<<","
    //     <<_Rank4TensorMateValues[0](1,1,2,2)<<","
    //     <<_Rank4TensorMateValues[0](1,1,3,3)<<","
    //     <<_Rank4TensorMateValues[0](1,1,2,3)<<","
    //     <<_Rank4TensorMateValues[0](1,1,1,3)<<","
    //     <<_Rank4TensorMateValues[0](1,1,1,2)<<endl
    //     <<_Rank4TensorMateValues[0](2,2,1,1)<<","
    //     <<_Rank4TensorMateValues[0](2,2,2,2)<<","
    //     <<_Rank4TensorMateValues[0](2,2,3,3)<<","
    //     <<_Rank4TensorMateValues[0](2,2,2,3)<<","
    //     <<_Rank4TensorMateValues[0](2,2,1,3)<<","
    //     <<_Rank4TensorMateValues[0](2,2,1,2)<<endl
    //     <<_Rank4TensorMateValues[0](3,3,1,1)<<","
    //     <<_Rank4TensorMateValues[0](3,3,2,2)<<","
    //     <<_Rank4TensorMateValues[0](3,3,3,3)<<","
    //     <<_Rank4TensorMateValues[0](3,3,2,3)<<","
    //     <<_Rank4TensorMateValues[0](3,3,1,3)<<","
    //     <<_Rank4TensorMateValues[0](3,3,1,2)<<endl
    //     <<_Rank4TensorMateValues[0](2,3,1,1)<<","
    //     <<_Rank4TensorMateValues[0](2,3,2,2)<<","
    //     <<_Rank4TensorMateValues[0](2,3,3,3)<<","
    //     <<_Rank4TensorMateValues[0](2,3,2,3)<<","
    //     <<_Rank4TensorMateValues[0](2,3,1,3)<<","
    //     <<_Rank4TensorMateValues[0](2,3,1,2)<<endl
    //     <<_Rank4TensorMateValues[0](1,3,1,1)<<","
    //     <<_Rank4TensorMateValues[0](1,3,2,2)<<","
    //     <<_Rank4TensorMateValues[0](1,3,3,3)<<","
    //     <<_Rank4TensorMateValues[0](1,3,2,3)<<","
    //     <<_Rank4TensorMateValues[0](1,3,1,3)<<","
    //     <<_Rank4TensorMateValues[0](1,3,1,2)<<endl
    //     <<_Rank4TensorMateValues[0](1,2,1,1)<<","
    //     <<_Rank4TensorMateValues[0](1,2,2,2)<<","
    //     <<_Rank4TensorMateValues[0](1,2,3,3)<<","
    //     <<_Rank4TensorMateValues[0](1,2,2,3)<<","
    //     <<_Rank4TensorMateValues[0](1,2,1,3)<<","
    //     <<_Rank4TensorMateValues[0](1,2,1,2)<<endl<<endl;

    cout<<"********************"<<endl;

        }
    }
    

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