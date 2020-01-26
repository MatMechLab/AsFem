#include "MaterialSystem/MaterialSystem.h"

double Sign(const double &x,const double &y){
    if(y>=0.0){
        return abs(x);
    }
    else{
        return -abs(x);
    }
}

double BrackPos(const double &x){
    return 0.5*(abs(x)+x);
}
double BrackPos(const double &x,const double &y){
    return 0.5*(x+Sign(x,y));
}
double BrackNeg(const double &x){
    return 0.5*(x-abs(x));
}
double BrackNeg(const double &x,const double &y){
    return 0.5*(x-Sign(x,y));
}

void MaterialSystem::ModifyElasticFracMate(const int &nDim,const double &t,const double &dt,
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

    if(InputParams.size()<5){
        cout<<"*** Error: for phase field fracutre,5 params are required!!***"<<endl;
        cout<<"***        you should give E,nu,Gc,L and viscosity !!!     ***"<<endl;
        Msg_AsFem_Exit();
    }
    // InputParams[0]--->E
    // InputParams[1]--->nu
    // InputParams[2]--->Gc
    // InputParams[3]--->L
    // InputParams[4]--->viscosity

    // so in matevalues, we store:
    // Gc, L, viscosity
    // Psi, Psi_pos, Psi_neg
    // Hist,dHistdc(should be 0 currently)
    
    _MateValues[0]=InputParams[2];// for Gc
    _MateValues[1]=InputParams[3];// for L
    _MateValues[2]=InputParams[4];// viscosity

    // here we use Rank4-1 to store the linear Cijkl, the final one should be rank4-0
    _Rank4TensorMateValues[1].SetToZeros();
    _Rank4TensorMateValues[1].SetFromEandNu(InputParams[0],InputParams[1]);
    const double lambda=_Rank4TensorMateValues[1](1,1,2,2);
    const double mu=_Rank4TensorMateValues[1](1,2,1,2);

    // calculate the elastic strain
    _Rank2TensorMateValues[2].SetToZeros();
    if(nDim==2){
        _Rank2TensorMateValues[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2TensorMateValues[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    // this is our strain
    _Rank2TensorMateValues[0]=(_Rank2TensorMateValues[2]+_Rank2TensorMateValues[2].Transpose())*0.5;
    
    RankTwoTensor Strain,StrainPos,StrainNeg,I;
    Strain=_Rank2TensorMateValues[0];
    double EpsTr=Strain.Trace();

    RankTwoTensor eigvec;
    double eigval[3];

    RankFourTensor ProjPos=_Rank2TensorMateValues[0].CalcPostiveProjTensor(eigval,eigvec);
    RankFourTensor I4Sym(RankFourTensor::InitIdentitySymmetric4);

    I.SetToIdentity();
    StrainPos=ProjPos*Strain;
    StrainNeg=Strain-StrainNeg;
    
    RankTwoTensor StressPos,StressNeg;

    StressPos=(lambda*BrackPos(EpsTr)*I+2*mu*StrainPos);
    StressNeg=lambda*BrackNeg(EpsTr)*I+2*mu*StressNeg;
    
    // Finally, we can get the stress as:
    double D=gpU[nDim];// for 2D: ux uy d
                       // for 3D: ux uy uz d
    if(D>0.9999) D=1.0;
    if(D<1.0e-4) D=0.0;
    _Rank2TensorMateValues[1]=StressPos*D*D+StressNeg;// Sigma
    _Rank2TensorMateValues[2]=StressPos*2.0*D;        // dSigma/dD

    _Rank4TensorMateValues[0]=D*D*(lambda*BrackPos(1,EpsTr)*I.CrossDot(I)+2*mu*ProjPos)
                             +lambda*BrackNeg(1,EpsTr)*I.CrossDot(I)+2*mu*(I4Sym-ProjPos);


    // use the fifth one to calculate the vonMises stress
    _Rank2TensorMateValues[4]=_Rank2TensorMateValues[1];
    // stress deviator tensor
    _Rank2TensorMateValues[4](1,1)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[4](2,2)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[4](3,3)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();

    // !!! Remember, the MateValues is already occupied by input params
    vector<double> epos(3), eneg(3);
    for(int i=0;i<3;++i){
        epos[i] = (abs(eigval[i]) + eigval[i])/2.0;//<lambda>+
        eneg[i] =-(abs(eigval[i]) - eigval[i])/2.0;//<lambda>-
    }
    double etr=0.0;
    for (int i=0;i<3;++i){
        etr+=eigval[i];// trace(lambda_i)
    }
    double etrpos= (abs(etr)+etr)/2.0;//<tr(eps)>+
    double etrneg=-(abs(etr)-etr)/2.0;
    double pval=0.0,nval=0.0;
    for (int i=0;i<3;++i){
        pval+=epos[i]*epos[i];//(<lambda>+)^{2}
        nval+=eneg[i]*eneg[i];//(<lamdda>-)^{2}
    }
    _MateValues[0]=InputParams[2];// for Gc
    _MateValues[1]=InputParams[3];// for L
    _MateValues[2]=InputParams[4];// viscosity
    _MateValues[3]= lambda*etrpos*etrpos/2.0+mu*pval;// Positve free energy
    _MateValues[4]=-lambda*etrneg*etrneg/2.0+mu*nval;// negative free energy
    _MateValues[5]=_MateValues[3]*D*D+_MateValues[4];// free energy

    // Important!!! update the history value
    if(_MateValues[3]>gpHistOld[0]){
        gpHist[0]=_MateValues[3];            //H
        gpHist[1]=0.0;      //dH/dD
        _Rank2TensorMateValues[3]=StressPos; // dH/dstrain
    }
    else{
        gpHist[0]=gpHistOld[0];
        gpHist[1]=0.0;
        _Rank2TensorMateValues[3]=0.0;  //dH/dstrain 
    }
    gpHist[1]=0.0; // dH/dD should always equal to 0
    // _Rank2TensorMateValues[3]=0.0;  //dH/dstrain 


    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    // sij is the deviator tensor of stress
    _MateValues[6]=sqrt(1.5*_Rank2TensorMateValues[4].DoubleDot(_Rank2TensorMateValues[4]));
    // hydrostatic stress
    _MateValues[7]=_Rank2TensorMateValues[1].Trace()/3.0;

    _Rank2TensorMateValues[4].SetToIdentity();//used in mechanics elements

    if(nDim==2){
        _MateValues[ 8]=_Rank2TensorMateValues[1](1,1);//sigxx
        _MateValues[ 9]=_Rank2TensorMateValues[1](2,2);//sigyy
        _MateValues[10]=_Rank2TensorMateValues[1](1,2);//sigxy

        _MateValues[11]=_Rank2TensorMateValues[0](1,1);//strxx
        _MateValues[12]=_Rank2TensorMateValues[0](2,2);//stryy
        _MateValues[13]=_Rank2TensorMateValues[0](1,2);//strxy
    }
    else if(nDim==3){
        _MateValues[ 8]=_Rank2TensorMateValues[1](1,1);//sigxx
        _MateValues[ 9]=_Rank2TensorMateValues[1](2,2);//sigyy
        _MateValues[10]=_Rank2TensorMateValues[1](3,3);//sigzz
        _MateValues[11]=_Rank2TensorMateValues[1](2,3);//sigyz
        _MateValues[12]=_Rank2TensorMateValues[1](1,3);//sigzx
        _MateValues[13]=_Rank2TensorMateValues[1](1,2);//sigxy

        _MateValues[14]=_Rank2TensorMateValues[0](1,1);//strxx
        _MateValues[15]=_Rank2TensorMateValues[0](2,2);//stryy
        _MateValues[16]=_Rank2TensorMateValues[0](3,3);//strzz
        _MateValues[17]=_Rank2TensorMateValues[0](2,3);//stryz
        _MateValues[18]=_Rank2TensorMateValues[0](1,3);//strzx
        _MateValues[19]=_Rank2TensorMateValues[0](1,2);//strxy
    }
}